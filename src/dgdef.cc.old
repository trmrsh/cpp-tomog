/*

!!begin
!!title  Compute gaussian default image
!!author T.R.Marsh
!!date   5 April 2001
!!root   dgdef
!!index  dgdef
!!descr  Computes gaussian default image
!!css   style.css
!!class  Doppler images
!!head1  ddef - computes gaussian default image

!!emph{dgdef} computes the gaussian default. This involves a fair few FFTs and can take a while.

!!head2 Invocation

dgdef input fwhm (gfwhm) output!!break

!!head2 Arguments

!!table
!!arg{ input  }{ input image (e.g. from !!emph{dinit}, '-' for standard input.}
!!arg{ fwhm   }{ image fwhm, pixels}
!!arg{ gfwhm  }{ gamma fwhm, pixels}
!!arg{ output }{ output file name '-' for standard output.}
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
#include "trm_dmap.h"

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("map",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("fwhm",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("gfwhm",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("output", Subs::Input::GLOBAL, Subs::Input::PROMPT);

    std::string infile;
    input.get_value("map",  infile, "input", "input file");
    Dmap map(infile);
    float fwhm;
    input.get_value("fwhm",   fwhm, 5.f, 0.f, FLT_MAX, "FWHM for blurring in Vx, Vy (pixels)");
    float gfwhm;
    if(map.ngamma() > 1)
      input.get_value("gfwhm", gfwhm, 5.f, 0.f, FLT_MAX, "FWHM for blurring in gamma (pixels)");
    std::string outfile;
    input.get_value("output", outfile, "output", "output file");

    float ibuf[map.size()], obuf[map.size()];
    map.get(ibuf);

    Tomog::gaussdef(ibuf,map.nwave(),map.ngamma(),map.nside(), fwhm, gfwhm, obuf);

    map.set(obuf);
    map.write(outfile);
  }

  catch(const Dmap::Dmap_Error& err){
    std::cerr << "Dmap_Error: " << err << std::endl;
    exit(EXIT_FAILURE);
  }
  catch(const std::bad_alloc&){
    std::cerr << "Memory allocation error" << std::endl;
    exit(EXIT_FAILURE);
  }
  catch(const std::string& err){
    std::cerr << "string exception: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}


