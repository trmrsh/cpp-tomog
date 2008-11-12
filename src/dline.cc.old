/*

!!begin
!!title  Adds a line to an image
!!author T.R.Marsh
!!created  30 June 2003
!!revised  30 June 2003
!!root   dline
!!index  dline
!!descr  Adds a line to an image
!!css   style.css
!!class  Fake
!!class  Doppler images
!!head1  dline - adds a line to an image

!!emph{dline} adds a line spot to an image. It is smeared in
a gaussian manner across all three dimensions. the user must
specify the start and end point as well as the FWHM.

!!head2 Invocation

dline map nwave x1 y1 z1 x2 y2 z2 height width output!!break

!!head2 Arguments

!!table
!!arg{ map    }{ input image.}
!!arg{ nwave  }{ the wavelength number of the image to add to.}
!!arg{ x1, y1, z1 }{ the three coordinates of the start of the line.}
!!arg{ x2, y2, z2 }{ the three coordinates of the end of the line.}
!!arg{ height }{ peak intensity.}
!!arg{ width  }{ FWHM of gaussian (km/s).}
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
#include "trm_vec3.h"
#include "trm_input.h"
#include "trm_tomog.h"
#include "trm_dmap.h"

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("map",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nwave",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("x1",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("y1",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("z1",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("x2",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("y2",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("z2",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("height", Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("width",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("output", Subs::Input::LOCAL,  Subs::Input::PROMPT);

    std::string infile;
    input.get_value("map",  infile, "input", "file to add spot to");
    Dmap map(infile);
    int nwave; 
    input.get_value("nwave", nwave,   1, 1, map.nwave(),  "which wavelength to operate on");
    nwave--;

    float x1;
    input.get_value("x1", x1, 0.f, -FLT_MAX, FLT_MAX, "X position of start of line (km/s)");
    float y1;
    input.get_value("y1", y1, 0.f, -FLT_MAX, FLT_MAX, "Y position of start of line (km/s)");
    float z1;
    input.get_value("z1", z1, 0.f, -FLT_MAX, FLT_MAX, "Z position of start of line (km/s)");
    float x2;
    input.get_value("x2", x2, 0.f, -FLT_MAX, FLT_MAX, "X position of end of line (km/s)");
    float y2;
    input.get_value("y2", y2, 0.f, -FLT_MAX, FLT_MAX, "Y position of end of line (km/s)");
    float z2;
    input.get_value("z2", z2, 0.f, -FLT_MAX, FLT_MAX, "Z position of end of line (km/s)");
    float height;
    input.get_value("height", height, 1.f, -FLT_MAX, FLT_MAX, "peak intensity of line");
    float width;
    input.get_value("width", width, 100.f, 0.0001f, FLT_MAX, "FWHM spread (km/s)");
    std::string outfile;
    input.get_value("output", outfile, "output", "output file");

    int nside = map.nside();
    float vpix   = map.vpix();

    Subs::Vec3 r1(x1,y1,z1);
    Subs::Vec3 r2(x2,y2,z2);
    Subs::Vec3 r, diff = r2 - r1;

    float efac = 1./Subs::sqr(width/Constants::EFAC)/2.;
    float cen  = float(nside-1)/2.;
    double lambda;

    for(int iz=0; iz<map.ngamma(); iz++){
      for(int iy=0; iy<nside; iy++){
	for(int ix=0; ix<nside; ix++){

	  // Set vector to pixel in question
	  r.set( vpix*(ix-cen), vpix*(iy-cen), map.gamma(iz));

	  // Find nearest point on the line connecting the two end-points to the pixel
	  lambda = dot(r - r1, diff)/diff.sqr();
	  lambda = std::max(0., std::min(1., lambda));
	  r -= (r1 + lambda*diff);
	  
	  // Add contribution according to how close the nearest point is.
	  map[nwave][iz][iy][ix] += height*exp(-efac*r.sqr());
	}
      }
    }

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


