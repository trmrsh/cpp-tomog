/*

!!begin
!!title  Generate a starting image
!!author T.R.Marsh
!!created  13 September 2000
!!revised  22 June 2003
!!root   dinit
!!index  dinit
!!descr  Generates a starting image
!!css   style.css
!!class  Fake
!!class  Doppler images
!!head1  dinit - generates a start image

!!emph{dinit} defines the basic structure of a Doppler image
file. Use !!ref{dspot.html}{dspot} to add in features.

!!head2 Invocation

dinit npix vpix wzero ngamma gamma (dgamma) output!!break

!!head2 Arguments

!!table
!!arg{ npix   }{ number of pixels on a side (same in x and y).}
!!arg{ vpix   }{ number of km/s/pixel.}
!!arg{ wzero  }{ vector of rest wavelengths.}
!!arg{ ngamma }{ number of gamma velocities.}
!!arg{ mgamma }{ mean gamma velocity (km/s).}
!!arg{ dgamma }{ spacing of gamma velocities (km/s) if ngamma > 1.}
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
    input.sign_in("nside",   Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("vpix",    Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("wzero",   Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("ngamma",  Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("mgamma",  Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("dgamma",  Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("output",  Subs::Input::LOCAL, Subs::Input::PROMPT);

    int nside;
    input.get_value("nside", nside, 100, 1, 2000, "number of pixels on a side of the images");
    float vpix;
    input.get_value("vpix",  vpix, 50.f, 0.001f, 1000.f, "pixel size (km/s)");
    std::vector<double> wzero;
    input.get_value("wzero", wzero, 5000., 1., 1.e6, 0, "wavelengths");
    int ngamma;
    input.get_value("ngamma",  ngamma, 100, 1, 2000, "number of systemic velocities");
    float mgamma;
    input.get_value("mgamma",  mgamma, 0.f, -FLT_MAX, FLT_MAX, "mean systemic velocity (km/s)");
    float dgamma;
    if(ngamma > 1)
      input.get_value("dgamma",  dgamma, 50.f, 0.001f, 1000.f, "step size in systemic velocity (km/s)");
    std::vector<float> gamma(ngamma);
    if(ngamma == 1){
      gamma[0] = mgamma;
    }else{
      for(int i=0; i<ngamma; i++)
	gamma[i] = mgamma+dgamma*(i-(ngamma-1)/2.);
    }
    std::string outfile;
    input.get_value("output", outfile, "output", "output file");

    // Create uninitialised image, set to zero, dump.
    Dmap dmap(nside,vpix,gamma,wzero);
    dmap = 0;
    dmap.write(outfile);
  }

  catch(const Dmap::Dmap_Error& err){
    std::cerr << "Dmap_Error: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  catch(const std::string& err){
    std::cerr << "string exception: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  catch(const std::bad_alloc&){
    std::cerr << "A memory allocation error has occurred!!" << std::endl;
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}


