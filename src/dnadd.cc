/*

!!begin
!!title  Add noise to Doppler images
!!author T.R.Marsh
!!created  09 July 2000
!!created  22 June 2003
!!root   dnadd
!!index  dnadd
!!descr  Adds noise to Doppler images
!!css   style.css
!!class  Fake
!!class  Noise
!!class  Doppler images
!!head1  dnadd - adds noise to Doppler images

!!emph{dnadd} adds noise to a Doppler image according to RMS values
defined by another. Its main use is to provide a comparison
to the correlated noise characteristic of true Doppler images.
The noise is added as gaussian white noise. See !!ref{dvar.html}{dvar}
for a routine that can generate the required RMS map.

!!head2 Invocation

dnadd input rms seed nout root!!break

!!head2 Arguments

!!table
!!arg{ map    }{ name of Doppler image to add noise to.}
!!arg{ rms    }{ name of Doppler image defining RMS values.}
!!arg{ seed   }{ random number seed integer.}
!!arg{ nout   }{ number of output files to generate.}
!!arg{ root   }{ root of output files. Thus if nout=100 the files will be
called names such as root001, root002 etc. If nout=1 then the output
is simply set = root.}
!!table

!!head2 Related commands

!!ref{tnadd.html}{tnadd}, !!ref{dvar.html}{dvar},
!!ref{tboot.html}{tboot}

!!end

 */

#include <cstdlib>
#include <climits>
#include <cmath>
#include <iostream>
#include <assert.h>
#include "trm_subs.h"
#include "trm_input.h"
#include "trm_tomog.h"
#include "trm_dmap.h"

void strint(const unsigned int num, const unsigned int nd, char* intstr);

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("map",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("rms",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("seed",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nout",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("root",    Subs::Input::LOCAL,  Subs::Input::PROMPT);

    std::string infile;
    input.get_value("map",  infile, "input", "input base file");
    std::string rmsfile;
    input.get_value("rms",  rmsfile, "rms", "rms file");
    Subs::INT4 seed;
    input.get_value("seed", seed,  57876, 1, INT_MAX, "seed integer for random number generator");
    seed = -seed;
    size_t nout;
    input.get_value("nout", nout,   size_t(1), size_t(1), size_t(1000000),  "the number of noisy files to be generated");
    std::string root;
    input.get_value("root",  root, infile, "root or base name for noisy files");

    // Read data in
    Dmap dmap(infile), rms(rmsfile), out;
    if(!match(dmap,rms))
      throw Tomog::Input_Error("Input and rms maps do not match");

    float sigma;
    size_t tnum = 10, ndg=1;
    while(tnum <= nout){
      tnum *= 10;
      ndg++;
    }

    for(size_t n=0; n<nout; n++){

      out = dmap;
    
      // add noise, dump to disc
      for(int nw=0; nw<out.nwave(); nw++){
	for(int ng=0; ng<out.ngamma(); ng++){
	  for(int iy=0; iy<out[nw][ng].nrow(); iy++){
	    for(int ix=0; ix<out[nw][ng].ncol(); ix++){
	      sigma = rms[nw][ng][iy][ix];
	      out[nw][ng][iy][ix] += sigma*Subs::gauss2(seed); 
	    }
	  }
	}
      }
      if(nout == 1){
	out.write(root);
      }else{
	out.write(root + "_" + Subs::str(int(n+1), ndg));
      }
    }
  }

  catch(const Dmap::Dmap_Error& err){
    std::cerr << "Dmap::Dmap_Error exception: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  catch(const std::string& err){
    std::cerr << "string exception: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}

  











