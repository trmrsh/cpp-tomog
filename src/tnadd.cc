/*

!!begin
!!title  Add noise to trailed spectra
!!author T.R.Marsh
!!date   13 September 2000
!!root   tnadd
!!index  tnadd
!!descr  Adds noise to trailed spectra
!!css   style.css
!!class  Fake
!!class  Noise
!!class  Trailed spectra
!!head1  tnadd - adds noise to trailed spectra
 
tnadd adds noise to a trailed spectrum and sets the
uncertainty array accordingly. It can be used
after !!ref{tgen.html}{tgen} to obtain a more
realistic noisy trail. tnadd can also simulate the effect
of poor slit losses by including overall multiplicative
factors to the spectra.

The noise on a pixel of flux = flux is generated according to 
sigma = sqrt(cont*(cont+flux))/ston, where cont is an effective 
continuum level, ston is the signal-to-noise in the ontinnum, 
and sigma is the RMS uncertainty.

!!head2 Invocation

tnadd trail cont ston slit seed nout root!!break

!!head2 Arguments

!!table
!!arg{ trail  }{ name of trailed spectrum.}
!!arg{ photon }{ photons/pixel/unit of flux.}
!!arg{ read   }{ readout noise/pixel.}
!!arg{ slit   }{ fractional RMS due to slit losses.}
!!arg{ seed   }{ random number seed integer.}
!!arg{ nout   }{ number of output files to generate.}
!!arg{ root   }{ root of output files. Thus if nout=10 the files will be
called names such as root1, root2 etc. If nout=1 then the output
is simply set = root.}
!!table

!!head2 Related commands

!!ref{dnadd.html}{tnadd}, !!ref{dvar.html}{dvar},
!!ref{tboot.html}{tboot}

!!end

 */

#include <climits>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <string>
#include <iostream>
#include "trm_subs.h"
#include "trm_input.h"
#include "trm_tomog.h"
#include "trm_trail.h"

void strint(const unsigned int num, const unsigned int nd, char* intstr);

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("trail",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("photon",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("read",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("slit",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("seed",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nout",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("root",    Subs::Input::LOCAL,  Subs::Input::PROMPT);

    std::string infile;
    input.get_value("trail",  infile, "input", "input base file");
    float photon;
    input.get_value("photon", photon, 1000.f, 0.e-5f, FLT_MAX,  "photons/pixels/flux unit");
    float read;
    input.get_value("read",   read, 10.f, 0.f, FLT_MAX,  "readout noise (photons/pixel)");
    float slit;
    input.get_value("slit", slit, 0.f, 0.f, 0.5f,  "RMS slit loss variability");
    Subs::INT4 seed;
    input.get_value("seed", seed,  57876, 1, INT_MAX, "seed integer for random number generator");
    seed = -seed;
    size_t nout;
    input.get_value("nout", nout,   size_t(1), size_t(1), size_t(1000000),  "the number of noisy files to be generated");
    std::string root;
    input.get_value("root",  root, infile, "root or base name for noisy files");

    Trail trail(infile);

    Trail tout;
    float sigma, mfac;
    size_t tnum = 10, ndg=1;
    while(tnum <= nout){
      tnum *= 10;
      ndg++;
    }

    for(unsigned int n=0; n<nout; n++){
      tout = trail;
    
      // add noise to it, dump to disc
      for(unsigned int i = 0; i < tout.nspec(); i++){
	mfac = std::max(0.,1.+slit*Subs::gauss2(seed));
	for(unsigned int j = 0; j < tout.npix(); j++){
	  tout.data()[i][j]  *= mfac;
	  sigma               = sqrt(read*read+photon*tout.data()[i][j])/photon;
	  tout.data()[i][j]  += sigma*Subs::gauss2(seed); 
	  tout.error()[i][j]  = sigma;
	}
      }
      if(nout == 1){
	tout.write(root);
      }else{
	tout.write(root + "_" + Subs::str(int(n+1), ndg));
      }
    }
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












