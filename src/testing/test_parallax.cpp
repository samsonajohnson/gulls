#include<cmath>
#include<vector>

#include "parallax.h"
#include "ephem.h"
#include "coords.h"
#include "constants.h"
#include "integerPowers.h"
#include "astroFns.h"

using namespace std;

int main()
{

  parallax pllx;
  parallax pllx2;
  parallax pllx3;
  coords c;
  vector<orbitalElements> orbit;
  
  double l=1.0*pi/180.0; //radians
  double b=-1.5*pi/180.0;
  double murel_l=5.0; //mas/yr
  double murel_b=-2.0; 
  double thetaE=1.0; //mas
  double piE=0.1;
  double tref=2460390.684641;

  //Setup the orbits 

  cout << scientific;

  orbit.resize(2);
  orbit[0].earthl2();
  orbit[1].lissajousz();

  pllx.reset();
  pllx.set_lb(l, b, 1);
  pllx.setup_reference_frame(tref,&orbit);
  pllx.set_orbit(&orbit);
  //pllx.setup_reference_frame();

  pllx.provide_murel_h_lb(murel_l, murel_b, piE, thetaE);

  pllx2.reset();
  pllx2.set_lb(l, b, 1);
  pllx2.setup_reference_frame(tref,&orbit);
  pllx2.set_orbit(&orbit);


  pllx2.provide_observables_NE(pllx.piEN, pllx.piEE, pllx.tE_r);

  cout << "tref" << " " << "a" << " " << "d" << " " << "l" << " " << "b" << " " << "lam" << " " << "bet" << " " << "mua_h" << " " << "mud_h" << " " << "mul_h" << " " << "mub_h" << " " << "mulam_h" << " " << "mubet_h" << " " << "murel_h" << " " << "mua_r" << " " << "mud_r" << " " << "mul_r" << " " << "mub_r" << " " << "mulam_r" << " " << "mubet_r" << " " << "murel_r" << " " << "vtilde_h" << " " << "vtilde_r" << " " << "v_ref" << " " << "vtilde_N_h" << " " << "vtilde_E_h" << " " << "vtilde_N_r" << " " << "vtilde_E_r" << " " << "v_ref_N" << " " << "v_ref_E" << " " << "piEN" << " " << "piEE" << " " << "piEll" << " " << "piErp" << " " << "phi_llN" << " " << "phi_pi" << " " << "piE" << " " << "thetaE" << " " << "rEtilde" << " " << "tE_h" << " " << "tE_r" << " " << endl;

  
  cout << tref << " " << pllx.a << " " << pllx.d << " " << pllx.l << " " << pllx.b << " " << pllx.lam << " " << pllx.bet << " " << pllx.mua_h << " " << pllx.mud_h << " " << pllx.mul_h << " " << pllx.mub_h << " " << pllx.mulam_h << " " << pllx.mubet_h << " " << pllx.murel_h << " " << pllx.mua_r << " " << pllx.mud_r << " " << pllx.mul_r << " " << pllx.mub_r << " " << pllx.mulam_r << " " << pllx.mubet_r << " " << pllx.murel_r << " " << pllx.vtilde_h << " " << pllx.vtilde_r << " " << pllx.v_ref << " " << pllx.vtilde_N_h << " " << pllx.vtilde_E_h << " " << pllx.vtilde_N_r << " " << pllx.vtilde_E_r << " " << pllx.v_ref_N << " " << pllx.v_ref_E << " " << pllx.piEN << " " << pllx.piEE  << " " << pllx.piEll << " " << pllx.piErp << " " << pllx.phi_llN << " " << pllx.phi_pi << " " << pllx.piE << " " << pllx.thetaE << " " << pllx.rEtilde << " " << pllx.tE_h << " " << pllx.tE_r << " " << endl;

  cout << tref << " " << pllx2.a << " " << pllx2.d << " " << pllx2.l << " " << pllx2.b << " " << pllx2.lam << " " << pllx2.bet << " " << pllx2.mua_h << " " << pllx2.mud_h << " " << pllx2.mul_h << " " << pllx2.mub_h << " " << pllx2.mulam_h << " " << pllx2.mubet_h << " " << pllx2.murel_h << " " << pllx2.mua_r << " " << pllx2.mud_r << " " << pllx2.mul_r << " " << pllx2.mub_r << " " << pllx2.mulam_r << " " << pllx2.mubet_r << " " << pllx2.murel_r << " " << pllx2.vtilde_h << " " << pllx2.vtilde_r << " " << pllx2.v_ref << " " << pllx2.vtilde_N_h << " " << pllx2.vtilde_E_h << " " << pllx2.vtilde_N_r << " " << pllx2.vtilde_E_r << " " << pllx2.v_ref_N << " " << pllx2.v_ref_E << " " << pllx2.piEN << " " << pllx2.piEE << " " << pllx2.piEll << " " << pllx2.piErp << " " << pllx2.phi_llN << " " << pllx2.phi_pi << " " << pllx2.piE << " " << pllx2.thetaE << " " << pllx2.rEtilde << " " << pllx2.tE_h << " " << pllx2.tE_r  << " " << endl;
  
  /*
  double l1, b1, l2, b2;
  double lam1, bet1, lam2, bet2;
  double a1, d1, a2, d2, a3, d3;
  double a, d;

  c.lb2ad(l,b,&a,&d);
  c.ad2lb(a,d,&l1,&b1);
  cout << "lb " << l*180.0/pi << " " << b*180.0/pi << endl;
  cout << "ad " << a*180.0/pi << " " << d*180.0/pi << endl;
  cout << "lb " << l1*180.0/pi << " " << b1*180.0/pi << endl;
  cout << endl << endl;

  c.ad2ecl(a,d,&lam1,&bet1);
  c.ecl2ad(lam1,bet1,&a2,&d2);
  cout << "ad " << a*180.0/pi << " " << d*180.0/pi << endl;
  cout << "\\s " << lam1*180.0/pi << " " << bet1*180.0/pi << endl;
  cout << "ad " << a2*180.0/pi << " " << d2*180.0/pi << endl;
  cout << endl << endl;

  eq2eclip(a*180./pi,d*180./pi,&lam1,&bet1);
  c.ecl2ad(lam1*pi/180.0,bet1*pi/180.0,&a2,&d2);
  //eclip2eq(lam1,bet1,&a2,&d2);
  cout << "ad " << a*180.0/pi << " " << d*180/pi << endl;
  cout << "\\s " << lam1 << " " << bet1 << endl;
  cout << "ad " << a2*180/pi << " " << d2*180/pi << endl;
  cout << endl << endl;


  c.mulb2ad(l,b,murel_l,murel_b,&a1,&d1);
  c.muad2lb(a,d,a1,d1,&l2,&b2);
  cout << "mulb " << murel_l << " " << murel_b << endl;
  cout << "muad " << a1 << " " << d1 << endl;
  cout << "mulb " << l2 << " " << b2 << endl;
  cout << endl << endl;
  
  c.muad2ecl(a,d,a1,d1,&lam2,&bet2);
  c.muecl2ad(a,d,lam2,bet2,&a3,&d3);
  cout << "muad " << a1 << " " << d1 << endl;
  cout << "mu\\s " << lam2 << " " << bet2 << endl;
  cout << "muad " << a3 << " " << d3 << endl;
  cout << endl << endl;
  */

  exit(0);

  for(tref=2460390.684641; tref<2460390.684641+366; tref+=1)
	{
	  pllx3.reset();
	  pllx3.set_lb(l, b, 1);
	  pllx3.setup_reference_frame(tref,&orbit);
	  pllx3.set_orbit(&orbit);

	  pllx3.provide_murel_h_lb(murel_l, murel_b, piE, thetaE);

	  cout << tref << " " << pllx3.a << " " << pllx3.d << " " << pllx3.l << " " << pllx3.b << " " << pllx3.lam << " " << pllx3.bet << " " << pllx3.mua_h << " " << pllx3.mud_h << " " << pllx3.mul_h << " " << pllx3.mub_h << " " << pllx3.mulam_h << " " << pllx3.mubet_h << " " << pllx3.murel_h << " " << pllx3.mua_r << " " << pllx3.mud_r << " " << pllx3.mul_r << " " << pllx3.mub_r << " " << pllx3.mulam_r << " " << pllx3.mubet_r << " " << pllx3.murel_r << " " << pllx3.vtilde_h << " " << pllx3.vtilde_r << " " << pllx3.v_ref << " " << pllx3.vtilde_N_h << " " << pllx3.vtilde_E_h << " " << pllx3.vtilde_N_r << " " << pllx3.vtilde_E_r << " " << pllx3.v_ref_N << " " << pllx3.v_ref_E << " " << pllx3.piEN << " " << pllx3.piEE << " " << pllx3.piEll << " " << pllx3.piErp << " " << pllx3.phi_llN << " " << pllx3.phi_pi << " " << pllx3.piE << " " << pllx3.thetaE << " " << pllx3.rEtilde << " " << pllx3.tE_h << " " << pllx3.tE_r  << " " << endl;

	}

  
}
