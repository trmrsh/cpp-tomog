/*

!!begin
!!title  Back-projects a trailed spectrum
!!author T.R.Marsh
!!created 13 September 2000
!!revised 22 June 2003
!!root   tback
!!index  tback
!!descr  Back-projects a trailed spectrum
!!css   style.css
!!class  Trailed spectra
!!class  Inversion
!!class  Doppler images
!!head1  tback - back-projects a trailed spectrum

!!emph{tback} computes the back-projection of a trailed spectrum.

!!head2 Invocation

tback trail tzero period nside vpix wzero gamma output!!break

!!head2 Arguments

!!table
!!arg{ trail  }{ trailed spectrum.}
!!arg{ tzero  }{ ephemeris zero-point.}
!!arg{ period }{ orbital period.}
!!arg{ nside  }{ number of pixels on a side in the output.}
!!arg{ vpix   }{ km/s/pixel.}
!!arg{ wzero  }{ rest wavelength.}
!!arg{ gamma  }{ systemic velocity (km/s).}
!!arg{ output }{ output image.}
!!table

!!end

*/

#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <iostream>
#include "trm_subs.h"
#include "trm_constants.h"
#include "trm_input.h"
#include "trm_tomog.h"
#include "trm_trail.h"
#include "trm_dmap.h"

int main(int argc, char *argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("trail",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("tzero",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("period",  Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("nside",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("vpix",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("wzero",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("gamma",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("output",  Subs::Input::LOCAL,  Subs::Input::PROMPT);

    std::string intrail;
    input.get_value("trail", intrail, "trail", "input trailed spectrum");
    Trail trail(intrail);
    double tzero;
    input.get_value("tzero",  tzero, 0.,  -DBL_MAX, DBL_MAX, "zero-crossing time");
    double period;
    input.get_value("period", period, 0.1, 1.e-6, DBL_MAX, "period");
    size_t nside;
    input.get_value("nside", nside, size_t(100), size_t(1), size_t(10000), "number of pixels along a side");
    float vpix;
    input.get_value("vpix", vpix, 50.f, 0.0001f, 10000.f, "number of km/s per pixel");
    double wzero;
    input.get_value("wzero", wzero, 5000., 0.0001, 1000000., "rest wavelength");
    float gamma;
    input.get_value("gamma", gamma, 0.f, -FLT_MAX, FLT_MAX, "systemic velocity (km/s)");
    std::string outfile;
    input.get_value("output", outfile, "map", "output Doppler map");

    // get trail info, compute cosines and sines
    int    nspec  = trail.nspec();
    int    npixd  = trail.npix();
    float  vpixd  = trail.vpix();
    double wzerod = trail.wzero();
    Subs::Array1D<double> cosp(nspec), sinp(nspec);

    cosp = cos(Constants::TWOPI*(trail.time()-tzero)/period);
    sinp = sin(Constants::TWOPI*(trail.time()-tzero)/period);

    // create map (single image)
    Dmap map(nside,vpix,gamma,wzero);

    int ilow, ns;
    float pvx, pvy, pvr, sum;
    float poff =  (Constants::C*1.e-3*(1.-wzerod/wzero)+gamma)/vpixd + (npixd-1)/2. + 0.5;

    for(size_t iy=0; iy < nside; iy++){
      pvy   = vpix*(iy-(nside-1)/2.)/vpixd;
      for(size_t ix=0; ix < nside; ix++){
	pvx   = vpix*(ix-(nside-1)/2.)/vpixd;
        for(ns=0, sum=0.; ns < nspec; ns++){
          pvr   = poff - pvx*cosp[ns] + pvy*sinp[ns];
          ilow  = int(floor(pvr));
	  if(ilow >= -1 && ilow < npixd){
	    if(ilow >= 0)      sum += (1+ilow-pvr)*trail.data()[ns][ilow];
	    if(ilow+1 < npixd) sum += (pvr-ilow)*trail.data()[ns][ilow+1];
	  }
	}
	map[0][0][iy][ix] = sum;
      }
    }

    map.write(outfile);

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
