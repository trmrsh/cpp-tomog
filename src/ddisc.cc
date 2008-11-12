/*

!!begin
!!title  Add a disc to an image
!!author T.R.Marsh
!!created  23 June 2003
!!revised  23 June 2003
!!root   dspot
!!index  dspot
!!descr  Add a disc to an image
!!css   style.css
!!class  Fake
!!class  Doppler images
!!head1  ddisc - adds a classic ring-like accretion disc to an image

!!emph{ddisc} adds an accretion disc represented as a ring. 

!!head2 Invocation

ddisc map nwave ngamma x0 y0 v1 v2 expon height gamma wgamma output!!break

!!head2 Arguments

!!table
!!arg{ map    }{ input image (e.g. from !!emph{dinit})}
!!arg{ nwave  }{ the wavelength of the image to add to.}
!!arg{ nwave  }{ the image to add to.}
!!arg{ x0     }{ x position of gaussian (km/s).}
!!arg{ y0     }{ y position of gaussian (km/s).}
!!arg{ v1     }{ the lowest velocity (km/s).}
!!arg{ v2     }{ the highest velocity (km/s).}
!!arg{ expon  }{ exponent of power law.}
!!arg{ height }{ height at lowest velocity.}
!!arg{ gamma  }{ mean systemic velocity (km/s).}
!!arg{ wgamma }{ FWHM of systemic velocity (km/s).}
!!arg{ output }{ output file name '-' for standard output.}
!!table

!!end

*/

#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <iostream>
#include "trm_constants.h"
#include "trm_subs.h"
#include "trm_input.h"
#include "trm_tomog.h"
#include "trm_dmap.h"
#include "trm_array2d.h"

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("map",    Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("nwave",  Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("x0",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("y0",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("v1",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("v2",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("expon",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("height", Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("gamma",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("wgamma", Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("output", Subs::Input::LOCAL,  Subs::Input::PROMPT);

    std::string infile;
    input.get_value("map",  infile, "input", "file to add spot to");
    Dmap map(infile);
    int nwave; 
    input.get_value("nwave", nwave,   1, 1, map.nwave(),  "which wavelength to operate on");
    nwave--;

    float x0;
    input.get_value("x0", x0, 0.f, -FLT_MAX, FLT_MAX, "X centre of gaussian spot (km/s)");
    float y0;
    input.get_value("y0", y0, 0.f, -FLT_MAX, FLT_MAX, "Y centre of gaussian spot (km/s)");
    float v1;
    input.get_value("v1", v1, 600.f, 0.f, FLT_MAX, "lowest velocity (outer disc, km/s)");
    float v2;
    input.get_value("v2", v2, 3000.f, v1, FLT_MAX, "highest velocity (inner disc, km/s)");
    float expon;
    input.get_value("expon", expon, -2.f, -30.f, 30.f, "exponent of intensity versus velocity");
    float height;
    input.get_value("height", height, 1.f, -FLT_MAX, FLT_MAX, "outer disc/low velocity intensity");
    float gamma;
    input.get_value("gamma", gamma, 0.f, -FLT_MAX, FLT_MAX, "systemic velocity of spot (km/s)");
    float wgamma;
    input.get_value("wgamma", wgamma, 100.f, 0.0001f, FLT_MAX, "FWHM of spot in systemic velocity terms (km/s)");
    std::string outfile;
    input.get_value("output", outfile, "output", "output file");

    int nside = map.nside();
    float vpix   = map.vpix();

    // Compute disc
    Subs::Array2D<float> disc(nside,nside);
    float cen  = float(nside-1)/2.;
    x0 = x0/vpix + cen;
    y0 = y0/vpix + cen;
    float v;
    for(int iy=0; iy<nside; iy++){
      for(int ix=0; ix<nside; ix++){
	v = vpix*sqrt(Subs::sqr(ix-x0) + Subs::sqr(iy-y0));

	if(v > v1 && v < v2){
	  disc[iy][ix] = height*pow(v/v1, expon);
	}else{
	  disc[iy][ix] = 0.;
	}

      }
    }

    // Add it in
    for(int ng=0; ng<map.ngamma(); ng++)
      map[nwave][ng] += float(exp(-Subs::sqr((map.gamma(ng)-gamma)/(wgamma/Constants::EFAC))/2.))*disc;

    map.write(outfile);
  }

  catch(const Dmap::Dmap_Error& err){
    std::cerr << "Dmap::Dmap_Error exception: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  catch(const std::string& err){
    std::cerr << "string exception: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  catch(const std::bad_alloc&){
    std::cerr << "Memory allocation error" << std::endl;
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}


