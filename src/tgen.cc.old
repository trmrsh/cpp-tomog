/*

!!begin
!!title  Generate a trailed spectrum
!!author T.R.Marsh
!!created 13 September 2000
!!revised 22 June 2003
!!root   tgen
!!index  tgen
!!descr  Compute/generate a trailed spectrum from a map
!!css   style.css
!!class  Trailed spectra
!!class  Doppler image
!!head1  tgen - compute a trailed spectrum from a map

!!emph{tgen} computes a trailed specturm equivalent to a map. It
does so by calling exactly the same routine as is used
by !!ref{dtmem.html}{dtmem} during the mem inversion. !!emph{tgen}
prompts for all the parameters needed to define a trail on a regular
set of orbital phases. It stores the phases as times 
effectively assuming tzero=0, period=1.

!!head2 Invocation

tgen map vpix wzero npix nspec phase1 phase2 ndiv fwhm file!!break

!!head2 Arguments

!!table
!!arg{ map     }{ name of a Doppler image.}
!!arg{ vpix    }{ km/s/pixel for trail.}
!!arg{ wzero   }{ central wavelength of trail. }
!!arg{ npix    }{ number of pixels/spectrum in the trail. }
!!arg{ nspec   }{ number of spectra in the trail. }
!!arg{ phase1  }{ first phase. }
!!arg{ phase2  }{ last phase. }
!!arg{ exposure}{ exposure time (in terms of orbital phase, used to simulate finite exposures) }
!!arg{ ndiv    }{ sub-division factor for improving projection accuracy. }
!!arg{ ntdiv   }{ number of points per spectrum to simulate finite exposure lengths }
!!arg{ fwhm    }{ fwhm (km/s) blurring. }
!!arg{ output  }{ output file name. }
!!table

!!end

 */

#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <iostream>
#include "trm_subs.h"
#include "trm_input.h"
#include "trm_array1d.h"
#include "trm_tomog.h"
#include "trm_dmap.h"
#include "trm_trail.h"

int main(int argc, char* argv[]){

  try{

        // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("map",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("vpix",    Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("wzero",   Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("npix",    Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("nspec",   Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("phase1",  Subs::Input::LOCAL,   Subs::Input::PROMPT);
    input.sign_in("phase2",  Subs::Input::LOCAL,   Subs::Input::PROMPT);
    input.sign_in("expose",  Subs::Input::LOCAL,   Subs::Input::PROMPT);
    input.sign_in("ndiv",    Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("ntdiv",   Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("fwhm",    Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("output",  Subs::Input::LOCAL,   Subs::Input::PROMPT);

    std::string inmap;
    input.get_value("map",     inmap,   "map",   "input Doppler map");
    float vpixd;
    input.get_value("vpix",    vpixd, 50.f,  0.0001f, 100000.f, "pixel size (km/s)");
    double wzerod;
    input.get_value("wzero",   wzerod, 5000., 0.0001, 1000000., "rest wavelength");
    int npixd;
    input.get_value("npix",    npixd, 100, 1, 100000, "number of pixels/spectrum");
    int nspec;
    input.get_value("nspec",   nspec, 100, 1, 100000, "number of spectra");
    double phase1;
    input.get_value("phase1",  phase1, 0.,  -DBL_MAX, DBL_MAX, "first phase");
    double phase2;
    input.get_value("phase2",  phase2, 1.,  -DBL_MAX, DBL_MAX, "final phase");
    float exposure;
    input.get_value("expose",  exposure, 0.01f, 0.f, 10.f, "exposure time (units of orbital period)");
    float fwhm;
    input.get_value("fwhm",    fwhm, 100.f, 0.0001f, 100000.f, "FWHM of local line profile (km/s)");
    int ndiv;
    input.get_value("ndiv",    ndiv, 1, 1, 200, "over-sampling factor for map/data computations");
    int ntdiv;
    input.get_value("ntdiv",   ntdiv, 1, 1, 200, "number of points/spectrum to simulate finite exposures");
    std::string outfile;
    input.get_value("output",  outfile, "trail", "output trail");

    Dmap map(inmap);

    // Create buffers for data and model. We work with these
    // rather than Trail and Dmap objects for compatibility
    // with the single-large-array nature of memsys.

    float mapbuf[map.size()];
    float datbuf[npixd*nspec];
    float vpix   = map.vpix();
    size_t nside = map.nside();

    // Transfer map data

    map.get(mapbuf);
    Subs::Array1D<float>  gamma = map.gamma();
    Subs::Array1D<double> wave  = map.wave();

    // Set times (effectively assumes ephemeris of tzero=0., period=1.)
    
    Subs::Array1D<double> time(nspec);
    Subs::Array1D<float> expose(nspec);
    for(int i=0; i<nspec; i++){
      if(nspec == 1){
	time[i] = (phase1+phase2)/2.;
      }else{
	time[i] = phase1 + (phase2-phase1)*i/(nspec-1.);
      }    
      expose[i] = exposure;
    }

    Tomog::op(mapbuf, wave, gamma, nside, vpix, fwhm, ndiv, ntdiv, npixd, 
		nspec, vpixd, wzerod, time, expose, 0., 1., datbuf);

    // Create and set trail
    Trail trail(npixd,nspec,vpixd,wzerod);
    trail.time()   = time;
    trail.expose() = expose;
    trail.set_data(datbuf);
    
    // set errors negative to indicate no noise
    for(unsigned int i = 0; i < trail.nspec(); i++){
      for(unsigned int j = 0; j < trail.npix(); j++){
	trail.error()[i][j]  = -1.;
      }
    }
    trail.write(outfile);
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











