#include "lightcurveGenerator.h"
#include "lightcurveFitter.h"
#include "astroFns.h"
#include "integerPowers.h"
#include "wittFSPL.h"
#include "singleLens.h"
#include "fisherinversion.h"

#define DEBUGVAR 1

void fisherMatrix(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses, int usefit)
{
  double t0, tE, u0;
  double tt,uu;
  vector<parallax> parlx;

  int idx,obsidx;
  
  if(DEBUGVAR) cout << "Calculating fisher matrix for event number " << Event->id << endl;

  int Nparams=3;
  int Nobsparams = 2*Paramfile->numobservatories;
  int Ntotparams = Nparams + Nobsparams;
  // 0 = t0  1 = tE  2 = u0  4 = F0  5 = fs
  Event->dF.resize(Ntotparams*Event->nepochs);
  double* fs = new double[Paramfile->numobservatories];
  double* F0 = new double[Paramfile->numobservatories];
  double u;

  vector<int> idxshift;
  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    idxshift.push_back(Event->nepochsvec[obsidx]);

  //first initialize the derivatives to zero so that we can accumulate the sum
  for(int param=0; param<Ntotparams; param++)
    {
      int pshift = param*Event->nepochs;
      for(idx=0; idx<Event->nepochs; idx++) Event->dF[pshift+idx]=0;
    }

  if(usefit>=0)
    {
      t0 = Event->PSPL[usefit].t0;
      tE = Event->PSPL[usefit].tE;
      u0 = Event->PSPL[usefit].umin;
      for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
	{
	  F0[obsidx] = Event->PSPL[usefit].Fl[obsidx] 
	    + Event->PSPL[usefit].Fu[obsidx];
	  fs[obsidx] = Event->PSPL[usefit].Fl[obsidx]/F0[obsidx];
	}
    }
  else
    {
      t0 = Event->t0;
      tE = Event->tE_r;
      u0 = Event->u0;
      for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
	{
	  fs[obsidx] = Event->fs[obsidx];
	  F0[obsidx] = 1;
	}
    }
  
  Event->deterror=0;

  double dmudu;
  double tau;
  int F0idx, fsidx;

  //Calculate the differential of the lightcurve
  for(idx=0;idx<Event->nepochs;idx++)
    {
      obsidx = Event->obsidx[idx];
      
      tau = (Event->epoch[idx] - t0) / tE;
      tt = tau;
      uu = u0;
            
      u = qAdd(tt,uu);

      dmudu = pspl_dmudu(u);
      double mu = pacAmp(u);
      
      // 0 = t0  1 = tE  2 = u0  4 = F0  5 = fs
      Event->dF[idx] = F0[obsidx] * fs[obsidx] * dmudu * tt/u * (-1.0/tE);
      Event->dF[idx+Event->nepochs] = F0[obsidx] * fs[obsidx] * dmudu * tt/u 
	* (-tau/tE);
      Event->dF[idx+2*Event->nepochs] = F0[obsidx] * fs[obsidx] * dmudu * uu/u;
      
      F0idx = (Nparams+obsidx*2)*Event->nepochs+idx;
      fsidx = F0idx + Event->nepochs;
      Event->dF[F0idx] = 1 + fs[obsidx]*(mu-1.0);
      Event->dF[fsidx] = F0[obsidx]*(mu-1.0);
      
    } //end for idx
  
  delete[] fs;
  delete[] F0;

}

