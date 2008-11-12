/*

!!begin
!!title  Add a gaussian to trailed spectra
!!author T.R.Marsh
!!date   25 September 2003
!!root   tgauss
!!index  tgauss
!!descr  Adds a gaussian to trailed spectra
!!css   style.css
!!class  Fake
!!class  Trailed spectra
!!head1  tgauss - adds a gaussian to trailed spectra

!!emph{tgauss} adds a gaussian of specifiable width, strength and 
radial velocity parameters to a trailed spectrum file. It can be
used multiple time to build up complex simulated trails. One can
also apply spectral broadening and smearing.

!!head2 Invocation

tgauss trail cont ston slit seed nout root!!break

!!head2 Arguments

!!table
!!arg{trail }{ name of input trailed spectrum.}
!!arg{height}{Height of gaussian}
!!arg{fwhm}{Full Width at Half Maximum in km/s if positive, Angstroms if negative}
!!arg{gamma}{Systemic velocity}
!!arg{kx}{Semi-amplitude in relation - kx*cos(2*pi*phi) + ky*sin(2*pi*phi)
where phi is the orbital phase. Sign ensure standard location in Doppler
images.}
!!arg{ky}{Semi-amplitude in relation - kx*cos(2*pi*phi) + ky*sin(2*pi*phi)
where phi is the orbital phase. Sign ensure standard location in Doppler
images.}
!!arg{t0}{Zero point of linear ephemeris}
!!arg{period}{Period of linear ephemeris}
!!arg{broad}{FWHM of spectral  broadening to apply, km/s}
!!arg{expose}{Length of exposures to account for smearing}
!!arg{ndiv}{Number of subdivisions for smearing}
!!arg{output }{ name of output trailed spectrum. Can be the same as the
input if wanted.}
!!table

!!end

 */

#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <string>
#include <iostream>
#include "trm_subs.h"
#include "trm_constants.h"
#include "trm_input.h"
#include "trm_tomog.h"
#include "trm_trail.h"

void strint(const unsigned int num, const unsigned int nd, char* intstr);

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("trail",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("height",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("fwhm",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("gamma",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("kx",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("ky",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("t0",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("period",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("broad",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("expose",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("ndiv",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("output",  Subs::Input::LOCAL,  Subs::Input::PROMPT);

    std::string infile;
    input.get_value("trail",  infile, "input", "input trailed spectrum");
    float height;
    input.get_value("height", height, 1.f, -FLT_MAX, FLT_MAX,  "gaussian peak height");
    float fwhm;
    input.get_value("fwhm",   fwhm, 100.f, -10000.f, 10000.f,  "FWHM of gaussian (+ve km/s, -ve Angstroms)");
    float gamma;
    input.get_value("gamma", gamma, 0.f, -10000.f, 10000.f,  "systemic velocity (km/s)");
    float kx;
    input.get_value("kx",    kx, 0.f, -10000.f, 10000.f,  "Kx value (km/s)");
    float ky;
    input.get_value("ky",    ky, 0.f, -10000.f, 10000.f,  "Ky value (km/s)");
    double t0;
    input.get_value("t0",    t0, 0., -DBL_MAX, DBL_MAX,  "zero point of ephemeris to set phases");
    double period;
    input.get_value("period", period, 1., 1.e-10, 1.e10,  "period of ephemeris to set phases");
    float broad;
    input.get_value("broad", broad, 100.f, 1.e-10f, 1.e5f,  "FWHM spectral  broadening to apply (km/s)");
    float expose;
    input.get_value("expose", expose, 0.1f, 0.f, 1.e5f,  "length of exposures for smearing");
    int ndiv;
    input.get_value("ndiv", ndiv, 1, 1, 1000,  "number of sub-divisions for smearing");
    std::string output;
    input.get_value("output",  output, infile, "output trailed spectrum");

    Trail trail(infile);

    // Apply spectral broadening
    if(fwhm < 0.) fwhm = - Constants::C/1000.*fwhm/trail.wzero();
    float fwhm_new = sqrt(Subs::sqr(fwhm) + Subs::sqr(broad));
    height *= (fwhm/fwhm_new);

    double sum[trail.npix()];
    const double sigma = fwhm_new/Constants::EFAC/trail.vpix();
    const double efac  = 1./Subs::sqr(sigma)/2.;
    double phi, v, wgt = 1.;
    for(unsigned int i = 0; i < trail.nspec(); i++){

      for(unsigned int j = 0; j < trail.npix(); j++)
	sum[j] = 0.;

      for(int k = 0; k < ndiv; k++){
	phi = Constants::TWOPI*(trail.time()[i]-t0)/period;
	if(ndiv > 1){
	  phi += expose*(k-(ndiv-1)/2.)/(ndiv-1);
	  if(k == 0 || k == ndiv-1) 
	    wgt = 0.5;
	  else
	    wgt = 1.;
	}
	v = (gamma - kx*cos(phi) + ky*sin(phi))/trail.vpix(); 
	for(unsigned int j = 0; j < trail.npix(); j++)
	  sum[j] += height*exp(-efac*Subs::sqr(j-(trail.npix()-1)/2.-v));
      }
      for(unsigned int j = 0; j < trail.npix(); j++)
	trail.data()[i][j] += sum[j]/std::max(1,ndiv-1);
    }

    trail.write(output);

  }

  catch(const Trail::Trail_Error& err){
    std::cerr << "Trail::Trail_Error exception: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  catch(const std::string& err){
    std::cerr << "string exception: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);

}












