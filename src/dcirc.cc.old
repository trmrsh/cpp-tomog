/*

!!begin
!!title  Convolution with a circle
!!author T.R.Marsh
!!created   13 September 2000
!!revised   21 June 2003
!!root   dcirc
!!index  dcirc
!!descr  Convolves an image with a circle
!!css   style.css
!!class  Arithematic
!!class  Doppler images
!!head1  dcirc - convoluion of a Doppler image with a circle.

dcirc convolves a Doppler image with a circle and produces a new
image which is the convolution. It is meant to allow calculations
to do with estimation of fluxes in apertures.

!!head2 Invocation

dcirc input r output!!break

!!head2 Arguments

!!table
!!arg{ input  }{ name of Doppler image file.}
!!arg{ r      }{ radius of circle to convolve image with.}
!!arg{ output }{ output Doppler image}
!!table

!!end

*/

#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <fstream>
#include "trm_subs.h"
#include "trm_constants.h"
#include "trm_input.h"
#include "trm_tomog.h"
#include "trm_dmap.h"

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("input",  Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("radius", Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("output", Subs::Input::GLOBAL, Subs::Input::PROMPT);

    std::string infile;
    input.get_value("input",  infile, "input", "input file");
    float radius;
    input.get_value("radius", radius, 5.f, 0.f, FLT_MAX, "convolution radius (km/s)");
    std::string outfile;
    input.get_value("output", outfile, "output", "output file");

    Dmap inmap(infile);
    Dmap output = inmap;

    int nsum;
    float sum, dysq;
    float vpix = inmap.vpix();
    int off    = int(radius/vpix) + 1;
    float rsq  = Subs::sqr(radius/vpix);
    for(int nw=0; nw<inmap.nwave(); nw++){
      for(int ng=0; ng<inmap.ngamma(); ng++){
	int nr   = inmap[nw][ng].nrow();
	int nc   = inmap[nw][ng].ncol();
	for(int iy=0; iy<nr; iy++){
	  for(int ix=0; ix<nc; ix++){
	    nsum = 0;
	    sum  = 0.;
	    for(int iy1=std::max(0,iy-off); iy1<std::min(nr,iy+off+1); iy1++){
	      dysq  = Subs::sqr(iy1 - iy);
	      for(int ix1=std::max(0,ix-off); ix1<std::min(nc,ix+off+1); ix1++){
		if(Subs::sqr(ix1-ix)+dysq < rsq){
		  sum += inmap[nw][ng][iy1][ix1];
		  nsum++;
		}
	      }
	    }
	    if(nsum){
	      output[nw][ng][iy][ix] = sum/nsum;
	    }else{
	      output[nw][ng][iy][ix] = 0.;
	    }
	  }
	}	 
      } 
    }
    output.write(outfile);
  }

  catch(const Dmap::Dmap_Error& err){
    std::cerr << "Dmap_Error: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}


