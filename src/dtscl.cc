/*

!!begin
!!title  Scales a map to minimise chi**2
!!author T.R.Marsh
!!date   13 September 2000
!!root   dtscl
!!index  dtscl
!!descr  Scales a map to minimise chi**2
!!css   style.css
!!class  Doppler images
!!class  Trailed spectra
!!head1  dtscl - scales a map to minimise chi**2

!!emph{dtscl} takes an image and a trailed spectrum and scales the
image to minimise the chi**2 relative to the data. This is
a simple overall scaling; at some point I will need to optimise
different lines separately.

!!head2 Invocation

dtscl map trail fwhm ndiv tzero period output!!break

!!head2 Arguments

!!table
!!arg{ map    }{ Doppler image.}
!!arg{ trail  }{ trailed spectrum.}
!!arg{ fwhm   }{ fwhm for blurring (km/s)}
!!arg{ ndiv   }{ sub-divison factor (>0)}
!!arg{ ntdiv  }{ time sub-divison factor (>0)}
!!arg{ tzero  }{ ephemeris zero-point.}
!!arg{ period }{ orbital period.}
!!arg{ output }{ output scaled image.}
!!table

!!end

*/

#include <cfloat>
#include <string>
#include "trm_subs.h"
#include "trm_input.h"
#include "trm_tomog.h"
#include "trm_dmap.h"
#include "trm_trail.h"

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("map",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("trail",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("fwhm",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("ndiv",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("ntdiv",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("tzero",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("period",  Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("output",  Subs::Input::LOCAL,  Subs::Input::PROMPT);

    std::string inmap;
    input.get_value("map",   inmap,   "map",   "input Doppler map");
    Dmap  dmap(inmap);
    std::string intrail;
    input.get_value("trail", intrail, "trail", "input trailed spectrum");
    Trail trail(intrail);
    float fwhm;
    input.get_value("fwhm", fwhm, 100.f, 0.0001f, 100000.f, "FWHM of local line profile (km/s)");
    int ndiv;
    input.get_value("ndiv", ndiv, 1, 1, 200, "over-sampling factor for map/data computations");
    int ntdiv;
    input.get_value("ntdiv", ntdiv, 1, 1, 200, "number of points per time point to simulate finite exposures");
    double tzero;
    input.get_value("tzero", tzero, 0.,  -DBL_MAX, DBL_MAX, "zero-crossing time");
    double period;
    input.get_value("period", period, 0.1, 1.e-6, DBL_MAX, "period");
    std::string outfile;
    input.get_value("output", outfile, "map", "output Doppler map");

    // Create and load buffers for data and model. 

    int npixd      = trail.npix();
    int nspec      = trail.nspec();
    float   vpix   = dmap.vpix();
    float   vpixd  = trail.vpix();
    double  wzerod = trail.wzero();

    float model[dmap.size()];
    dmap.get(model);
    Subs::Array1D<float>  gamma  = dmap.gamma();
    Subs::Array1D<double> wave   = dmap.wave();
    Subs::Array1D<double> time   = trail.time();
    Subs::Array1D<float>  expose = trail.expose();

    size_t ndat       = trail.size();
    float data[ndat], errors[ndat], calc[ndat];
    size_t nside = dmap.nside();
 
    Tomog::op(model, wave, gamma, nside, vpix, fwhm, ndiv, ntdiv, npixd, 
		nspec, vpixd, wzerod, time, expose, tzero, period, calc);

    trail.get_data(data);
    trail.get_error(errors);

    float sum1 = 0., sum2 = 0., chi1=0.;
    for(size_t i=0; i<ndat; i++){
      sum1 += data[i]*calc[i]/Subs::sqr(errors[i]);
      sum2 += Subs::sqr(calc[i]/errors[i]);
      chi1 += Subs::sqr((data[i]-calc[i])/errors[i]);
    }
    float scale = sum1/sum2, chi2=0.;
    for(size_t i=0; i<ndat; i++)
      chi2 += Subs::sqr((data[i]-scale*calc[i])/errors[i]);
    dmap *= scale;
    dmap.write(outfile);
    std::cerr << "Map scaled by factor = " << scale << std::endl;
    std::cerr << "Chi**2/N changed from " << chi1/ndat << " to " << chi2/ndat << std::endl;

  }

  catch(const Dmap::Dmap_Error& err){
    std::cerr << "Dmap::Dmap_Error exception: " << err << std::endl;
    exit(EXIT_FAILURE);
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


