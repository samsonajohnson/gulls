#include "lightcurveGenerator.h"
#include "backupGenerator.h"
#include "astroFns.h"
#include "VBBinaryLensingLibrary.h"
#include "fs.h"
#include "src_cld.h"
#include "lens_binary.h"
#include "cd.h"

#include<fstream>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>

#define DEBUGVAR 1

extern "C"
{
  void magfunc_(double *m1, double *a, double *xsCenter,  double *ysCenter, double *rs, double *Gamma, double *amp, double *eps, int *errflag);
}

void fisherMatrix(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses)
{
  double m1, a;
  double xsCenter, ysCenter, rs, Gamma;
  double amp, eps=1.0e-3;
  double alpha, cosa, sina, xcom;
  double t0, tE, u0;
  double tt,uu;
  double logq=0, logs=0, logtE=0;
  double pllx=Paramfile->pllxMultiplyer;
  vector<parallax> parlx;
  double lcgen=Paramfile->LC_GEN;
  int idx,obsidx;
  int useVBB=1;
  VBBinaryLensing VBBL; 
  int filter;

  int errflag;

  double piEN=0, piEE=0, piE;

  if(DEBUGVAR) cout << "Calculating fisher matrix for event number " << Event->id << endl;

  int Nparams=9;
  if(pllx==0) Nparams=7;
  int Nobsparams = 2*Paramfile->numobservatories;
  int Ntotparams = Nparams + Nobsparams;
  // 0 = t0  1 = tE  2 = u0  3 = alpha  4 = s  5 = q  6 = rs  7 = piEN 8 = piEE 9 = F0  10 = fs
  double* step = new double[Ntotparams];
  //double* dA = new double[2*MAX_NUM_EPOCH];
  double* fs = new double[Paramfile->numobservatories];
  int dshift;
  double dmult;
  double z1,z2;

  int pshift;
  if(DEBUGVAR) cout << "LC_GEN: " << lcgen << endl;

  src_cld src;
  lens_binary lens;
  finiteSource finsrc;
  src_cld altsrc;
  lens_binary altlens;
  finiteSource altfs;

  step[0] = 1e-6; step[1] = 1e-6; step[2] = 1e-6; step[3] = 1e-6;
  step[4] = 1e-6; step[5] = 1e-6; step[6] = 1e-6; 
  if(pllx!=0)
    {
      step[7] = 1e-6;
      step[8] = 1e-6;
    }
  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      step[Nparams+2*obsidx] = 1e-6;
      step[Nparams+2*obsidx+1] = 1e-6;
      parlx.push_back(Event->pllx[obsidx]);
    }

  Event->dF.resize(Ntotparams*Event->nepochs);
  //first initialize the derivatives to zero so that we can accumulate the sum
  for(int param=0; param<Ntotparams; param++)
    {
      int pshift = param*Event->nepochs;
      for(idx=0; idx<Event->nepochs; idx++) Event->dF[pshift+idx]=0;
    }

  //Calculate derivatives numerically for the non-linear parameters
  for(int param=0;param<Nparams;param++)
    {
      pshift = param*Event->nepochs;

      if(DEBUGVAR) cout << "Parameter: " << param << endl;
      for(int delta=-1;delta<=+1;delta+=2)
	{
	  
	  dshift = Event->nepochs*(delta+1)/2;
	  dmult = 0.5*delta/step[param];
	  //work out the event parameters in the fortran parametrization

	  //mass of the first lens m1+m2=1
	  logq = log10(Event->params[QQ]) + (param==5?delta*step[5]:0);
	  logs = log10(Event->params[SS]) + (param==4?delta*step[4]:0);
	  logtE = log10(Event->tE_r) + (param==1?delta*step[1]:0);

	  m1 = 1.0/(1 + pow(10,logq)); 
	  a = pow(10,logs);	//separation
	  z1 = -m1*a;
	  z2 = (1-m1)*a;
	  rs = pow(10,log10(Event->rs) + (param==6?delta*step[6]:0));	        
	  alpha = Event->alpha*TO_RAD + (param==3?delta*step[3]:0); 
	  Gamma = Event->gamma;	/* limb-darkening profile */
	  
	  if(pllx)
	    {
	      piEN = Event->piEN + (param==7?delta*step[7]:0);
	      piEE = Event->piEE + (param==8?delta*step[8]:0);
	      piE = qAdd(piEN,piEE);
	    }

	  //only change the finite source parameters when needed to aid efficiency
	  if(param==5 || param==4 || param==6)
	    {
	      altsrc = src_cld(rs,0);
	      altlens = lens_binary(m1, z1, 1-m1, z2);
	      altfs = finiteSource(&altsrc,&altlens,5.0e-4,1.0e-3);
	      altfs.memoryForRays(500);
	      altfs.saveRays();
	      altfs.reset();
	    }
	  else if(param==0 && delta==-1)
	    {
	      //This is the nominal values
	      src = src_cld(rs,0);
	      lens = lens_binary(m1, z1, 1-m1, z2);
	      finsrc = finiteSource(&src,&lens,5.0e-4,1.0e-3);
	      finsrc.memoryForRays(500);
	      finsrc.saveRays();
	      finsrc.reset();
	    }

	  cosa = cos(alpha); sina = sin(alpha);
	  /*xcom = a*(1-2*m1);*/
	  xcom = -m1*a;  /* Origin is primary lens*/

	  t0 = Event->t0 + (param==0?delta*step[0]:0);
	  tE = pow(10,logtE);
	  u0 = Event->u0 + (param==2?delta*step[2]:0);
	  //if(abs(Event->u0)<1.0e-5) u0 = Event->u0 + (param==2?delta*step[2]:0);
	  //else u0 = pow(10,log10(Event->u0) + (param==2?delta*step[2]:0));

	  //only adjust parallax parameters if necessary
	  if(pllx && (param==7 || param==8 || param == 1))
	    {
	      for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
		{
		  parlx[obsidx].set_tE_r(tE);
		  parlx[obsidx].set_piENE(piEN,piEE);
		  parlx[obsidx].fit_reinit();
		}
	    }


	  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
	    {
	      //fs[obsidx] = pow(10, log10(Event->fs[obsidx]) + (param==Nparams+2*obsidx+1?delta*step[Nparams+2*obsidx+1]:0));
	      fs[obsidx] = Event->fs[obsidx];
	    }

	  Event->deterror=0;
	  errflag=0;

	  //Calculate the lightcurve
	  for(idx=0;idx<Event->nepochs;idx++)
	    {
	      obsidx = Event->obsidx[idx];
	      
	      tt = (Event->epoch[idx] - t0) / tE;
	      uu = u0;

	      if(pllx)
		{
		  //tt += parlx[obsidx].tshift(idx-idxshift[obsidx]);
		  //uu += parlx[obsidx].ushift(idx-idxshift[obsidx]);
		  tt += parlx[obsidx].tshift(Event->jdepoch[idx]);
		  uu += parlx[obsidx].ushift(Event->jdepoch[idx]);
		}

	      xsCenter = tt*cosa - uu*sina + xcom;
	      ysCenter = tt*sina + uu*cosa;
              if(lcgen==0) {
                   magfunc_(&m1, &a, &xsCenter, &ysCenter, &rs, &Gamma, &amp, &eps,
                       &errflag);
                  // cout <<"old lcgen"<<endl;
               } else if (lcgen==1) {
                  VBBL.a1 = 0.36; // Linear limb darkening coefficient. 
                  amp = VBBL.BinaryMag2(m1, a, xsCenter, ysCenter,rs);
                  //cout <<"VBBL lcgen"<<endl;
              //xsCenter-=xcom; amp = pacAmp(qAdd(xsCenter,ysCenter));//for testing          
              } //end for lc_gen


              //if(param==5 || param==4 || param==6)
	      //	amp = altfs.fsmag(cd(xsCenter,ysCenter), errflag);
	      //else
	      //	amp = finsrc.fsmag(cd(xsCenter,ysCenter), errflag);	      
 
	      if( errflag != 0) 
		{
		  cerr << "\nerror caught from magfunc_  errval: " 
		       << Event->lcerror << " (fisher)" << endl;
		  errorHandler(errflag);
		  Event->deterror=errflag;
		  break;
		}

	      filter = World[obsidx].filter;

	      //store the results - remember F0=1 by definition
	      //	      dA[dshift+idx] = 1 + fs[obsidx]*(amp-1);
	      //dF = dA/dp = fs[obsidx]*amp
	      //dF[pshift+idx] += dmult*(1 + fs[obsidx]*(amp-1));
	      Event->dF[pshift+idx] += dmult*fs[obsidx]*amp;
  
	    } //end for idx
	} //end for delta

      //reset parallax parameters back to nominal values if necessary
      if(pllx && (param==7 || param==8 || param == 1))
	{
	  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
	    {
	      parlx[obsidx].set_tE_r(Event->tE_r);
	      parlx[obsidx].set_piENE(Event->piEN,Event->piEE);
	      parlx[obsidx].fit_reinit();
	    }
	}

      //compute dF/dparam

      //for(idx=0;idx<Event->nepochs;idx++)
      //	{
      //	  dF[pshift+idx] = 0.5 * (dA[MAX_NUM_EPOCH + idx] - dA[idx]) 
      //	    / (step[param]);
      //	}
    } //end for parameter

  //For the linear flux parameters we can work analytically
  //Do this numerically for the non-linear parameters
  int F0idx, fsidx;

  //Now calculate dF/dF0 and dF/dfs
  for(idx=0;idx<Event->nepochs;idx++)
    {
      obsidx = Event->obsidx[idx];
      F0idx = (Nparams+obsidx*2)*Event->nepochs+idx;
      fsidx = F0idx + Event->nepochs;
      //Because Atrue=fs*mu+1-fs
      Event->dF[F0idx]=Event->Atrue[idx];
      Event->dF[fsidx]=(Event->Atrue[idx]-1)/fs[obsidx]; 
    }     

  delete[] step;
  delete[] fs;

}

