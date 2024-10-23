#include<fstream>
#include<iostream>
#include<string>
#include<cmath>
#include<vector>
#include "integerPowers.h"
#include "constants.h"

using namespace std;

int main(int argc, char* argv[])
{
  if(argc!=6&&argc!=5)
    {
      cerr << "Usage: ./generatePSF <sampling (arcsec)> <nsamples> <fwhm> {<beta>} <output>" << endl;
      exit(1);
    }

  double sampling=atof(argv[1]);
  int nsamples=atoi(argv[2]);
  double fwhm=atof(argv[3]);
  double beta;
  string outname=string(argv[argc-1]);

  if(argc==6) beta=atof(argv[argc-2]);
  else beta=4;

  ofstream out(outname.c_str());
  if(!out)
    {
      cerr << "Could not open output file: " << outname << endl;
      exit(1);
    }

  double I;
  double r;
  double alpha = fwhm/(2.0*sqrt(pow(2.0,1.0/beta)-1));
  double norm = (beta-1.0)/(pi*sqr(alpha));

  cout << "alpha = " << alpha << "\n";
  cout << "beta = " << beta << "\n";
  cout << "norm = " << norm << endl;

  //for each pixel
  for(int i=-nsamples/2;i<=nsamples/2;i++)
    {
      for(int j=-nsamples/2;j<=nsamples/2;j++)
	{	  
	  r = sampling*qAdd(double(i),double(j))/alpha;
	  I = norm*pow(1.0 + r*r,-beta);

	  out << I << " ";
	    
	  
	} //end for j
      out << endl;
    } //end for i
}
