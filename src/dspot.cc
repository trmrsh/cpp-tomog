/*

!!begin
!!title  Add a gaussian spot to an image
!!author T.R.Marsh
!!created  5 April 2001
!!revised  22 June 2003
!!root   dspot
!!index  dspot
!!descr  Add a gaussian spot to an image
!!css   style.css
!!class  Fake
!!class  Doppler images
!!head1  dspot - adds a gaussian spot to image

!!emph{dspot} adds a gaussian spot to an image. It is gaussian
across the Vx, Vy coordinates of the images and along the systemic
velocity dorection.

!!head2 Invocation

dspot map nwave x0 y0 height width gamma wgamma output!!break

!!head2 Arguments

!!table
!!arg{ map    }{ input image (e.g. from !!emph{dinit}, '-' for standard input.}
!!arg{ nwave  }{ the image to add to.}
!!arg{ x0     }{ x position of gaussian (km/s).}
!!arg{ y0     }{ y position of gaussian (km/s).}
!!arg{ height }{ height of gaussian.}
!!arg{ width  }{ FWHM of gaussian (km/s).}
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
#include "trm_buffer2d.h"

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("map",    Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("nwave",  Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("x0",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("y0",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("height", Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("width",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
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
    float height;
    input.get_value("height", height, 1.f, -FLT_MAX, FLT_MAX, "height of gaussian spot");
    float width;
    input.get_value("width", width, 100.f, 0.0001f, FLT_MAX, "FWHM of spot in X and Y (km/s)");
    float gamma;
    input.get_value("gamma", gamma, 0.f, -FLT_MAX, FLT_MAX, "systemic velocity of spot (km/s)");
    float wgamma;
    input.get_value("wgamma", wgamma, 100.f, 0.0001f, FLT_MAX, "FWHM of spot in systemic velocity terms (km/s)");
    std::string outfile;
    input.get_value("output", outfile, "output", "output file");

    int nside = map.nside();
    float vpix   = map.vpix();

    Subs::Array2D<float> spot(nside,nside);
    float efac = 1./Subs::sqr(width/vpix/Constants::EFAC)/2.;
    float cen  = float(nside-1)/2.;
    x0 = x0/vpix + cen;
    y0 = y0/vpix + cen;
    for(int iy=0; iy<nside; iy++)
      for(int ix=0; ix<nside; ix++)
	spot[iy][ix] = height*exp(-efac*(Subs::sqr(ix-x0) + Subs::sqr(iy-y0)));

    for(int ng=0; ng<map.ngamma(); ng++)
      map[nwave][ng] += float(exp(-Subs::sqr((map.gamma(ng)-gamma)/(wgamma/Constants::EFAC))/2.))*spot;

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


