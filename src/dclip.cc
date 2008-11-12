/*

!!begin
!!title  Clip Doppler images
!!author T.R.Marsh
!!date   03 September 2000
!!root   dclip
!!index  dclip
!!descr  Thresholds Doppler images
!!css   style.css
!!class  Arithematic
!!class  Doppler images
!!head1  dclip - thresholds a Doppler image.

dclip sets all values above and below user-defined thresholds
to those values. It can be used for instance to ensure that there
are no negative values in an image.

!!head2 Invocation

dclip input lo hi output!!break

!!head2 Arguments

!!table
!!arg{ input }{ name of Doppler image file.}
!!arg{ lo    }{ lower threshold.}
!!arg{ hi    }{ upper threshold.}
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
#include "trm_input.h"
#include "trm_tomog.h"
#include "trm_dmap.h"

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("input",  Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("low",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("high",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("output", Subs::Input::LOCAL,  Subs::Input::PROMPT);

    std::string infile;
    input.get_value("input",  infile, "input", "input file");
    float low;
    input.get_value("low",  low,  0.f, -FLT_MAX, FLT_MAX, "lower clip limit");
    float high;
    input.get_value("high", high, std::max(low, 1000.f), low, FLT_MAX, "upper clip limit");
    std::string outfile;
    input.get_value("output", outfile, "output", "output file");

    Dmap map(infile);

    int nlo = 0, nhi = 0;
    for(int nw=0; nw<map.nwave(); nw++){
      for(int ng=0; ng<map.ngamma(); ng++){
	for(int iy=0; iy<map[nw][ng].nrow(); iy++){
	  for(int ix=0; ix<map[nw][ng].nrow(); ix++){
	    if(map[nw][ng][iy][ix] < low){
	      nlo++;
	      map[nw][ng][iy][ix] = low;
	    }else if(map[nw][ng][iy][ix] > high){
	      nhi++;
	      map[nw][ng][iy][ix] = high;
	    }
	  }
	}
      }
    }
    std::cout << nlo << " pixels were set = " << low  << std::endl;
    std::cout << nhi << " pixels were set = " << high << std::endl;
    map.write(outfile);
  }

  catch(const Dmap::Dmap_Error& err){
    std::cerr << "Dmap_Error: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  catch(const std::string& err){
    std::cerr << "Unrecognised std::string exception: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}


