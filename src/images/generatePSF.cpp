#include<fstream>
#include<iostream>
#include<string>
#include<cmath>
#include<vector>

#include<gsl/gsl_sf_bessel.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_sf_trig.h>

#include "integerPowers.h"
#include "constants.h"

using namespace std;

/*struct airyparams{
  double r;
  double diameter;
  };*/

struct asparams{
  double x,y; //position in focal plane in arcsec
  double diameter; //telescope diameter in m
  double obstruction; //obstruction diameter in m
  double W; //width of spider vanes
  double f; //telescope focal length
  double rot; //vane rotation
};

double airy(double x, void* params)
{
  asparams* a = (asparams*) params;
  double r=qAdd(a->x,a->y);
  double diameter=a->diameter;
  double lambda=x;

  double theta = pi*diameter*r*arcsec2rad/lambda;
  
  if(r<1.0e-10) return 1.0;
  else return sqr(2.0*gsl_sf_bessel_J1(theta)/theta);
}

double obsairy(double x, void* params)
{
  asparams* a = (asparams*) params;
  double r=qAdd(a->x,a->y);
  double obs=a->obstruction;
  double diam=a->diameter;
  double od2=sqr(obs/diam);
  double lambda=x;

  double thetad = pi*diam*r*arcsec2rad/lambda;
  double thetao = pi*obs*r*arcsec2rad/lambda;
  
  if(r<1.0e-10) return 1.0;
  else return 4*sqr(gsl_sf_bessel_J1(thetad)/thetad 
			       - od2*gsl_sf_bessel_J1(thetao)/thetao);
}

double spider(double x, void* params)
{
  asparams* p = (asparams*) params;
  
  double rot = p->rot;
  double crot=cos(rot);
  double srot=sin(rot);

  double xx = p->x*crot - p->y*srot;
  double yy = p->x*srot + p->y*crot;
  double W = p->W; //vane width, meters
  double diameter = p->diameter;
  double f = p->f; //focal length, meters

  double lambda=x;

  xx *= f*arcsec2rad/lambda; //convert angular position to position in terms of wavelengths.
  yy *= f*arcsec2rad/lambda; 

  //sinc is normalized sinc i.e. sinc(x)=sin(pi x)/(pi x)

  double bwx = pi*W/f; //pi already accounted for
  double bwy = -bwx*yy;
  bwx *= xx;
  double bhx = pi*diameter/f; //pi already accounted for
  double bhy = bhx*yy;
  bhx *= xx;

  return 0.5*(sqr(gsl_sf_sinc(bwx)*gsl_sf_sinc(bhy)) + sqr(gsl_sf_sinc(bwy)*gsl_sf_sinc(bhx)));
  
}

double airyspider(double x, void* params)
{
  asparams* p = (asparams*) params;

  double d = p->diameter;
  double W = p->W;
  double f = p->f;
  
  double lambda = x;

  double Is = 2*pi*d*W*W/(f*sqr(lambda));
  double Id = sqr(0.25*pi*d*d/lambda/f)-Is;

  /*double air = airy(x,params);
  double oair = obsairy(x,params);
  double spid = spider(x,params);

  cout << od << " " << Is << " " << Id << endl;
  cout << air << " " << od*od*oair << " " << spid << endl;
  */
  

  double val = Id*obsairy(x,params) + Is*spider(x,params);

  //cout << val << endl;

  return val;
}

int main(int argc, char* argv[])
{
  if(argc!=10&&argc!=11)
    {
      cerr << "Usage: ./generatePSF <lambda1 (um)> <lambda2 (um)> <diameter (m)> <obstruction (m)> <vane width (m)> <focal length (m)> <sampling (arcsec)> <nsamples> {<vane rotation=45 deg>} <output>" << endl;
      exit(1);
    }

  double lambda1=atof(argv[1])*1e-6;
  double lambda2=atof(argv[2])*1e-6;
  double diameter=atof(argv[3]);
  double obstruction=atof(argv[4]);
  double width=atof(argv[5]);
  double focal=atof(argv[6]);
  double sampling=atof(argv[7]);
  int nsamples=atof(argv[8]);
  double rot=45*pi/180;;
  if(argc==11) rot = atof(argv[argc-2])*pi/180;
  string outname=string(argv[argc-1]);

  asparams aparams;
  aparams.diameter=diameter;
  aparams.obstruction=obstruction;
  aparams.W=width;
  aparams.f=focal;
  aparams.rot=rot;

  ofstream out(outname.c_str());
  if(!out)
    {
      cerr << "Could not open output file: " << outname << endl;
      exit(1);
    }

  gsl_integration_workspace* W =  gsl_integration_workspace_alloc(1000);
  
  gsl_function F;
  F.function = &airyspider;
  F.params = &aparams;

  double I,err;
  double dlambda=lambda2-lambda1;

  //for each pixel
  for(int i=-nsamples/2;i<=nsamples/2;i++)
    {
      //#pragma omp parallel for
      for(int j=-nsamples/2;j<=nsamples/2;j++)
	{
	  aparams.x = sampling*double(i);
	  aparams.y = sampling*double(j);
	  //cout << i << " " << j << endl;

	  //if(i==0&&j==0) out << 1.0 << " ";
	  //else
	  //{
	      gsl_integration_qag(&F,lambda1,lambda2,1e-10,1.0e-5,1000,GSL_INTEG_GAUSS61,W,&I,&err);

	      out << I/dlambda << " ";
	      //}
	  
	} //end for j
      out << endl;
    } //end for i

  gsl_integration_workspace_free(W);
}
