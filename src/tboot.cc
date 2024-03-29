/*

!!begin
!!title  Bootstrap trail generation
!!author T.R.Marsh
!!created 13 September 2000
!!revised 22 June 2003
!!root   tboot
!!index  tboot
!!descr  Bootstrap generation of trailed spectra
!!css   style.css
!!class  Fake
!!class  Noise
!!head1  tboot - generate trailed spectra using the bootstrap method.

!!emph{tboot} generates multiple trailed spectra by bootstrap resampling.
It does this by modifying the uncertainties on the data.

!!head2 Invocation

tboot trail seed nout root!!break

!!head2 Arguments

!!table
!!arg{ trail  }{ name of initial trailed spectrum.}
!!arg{ seed   }{ seed integer.}
!!arg{ nout   }{ number of trails to generate.}
!!arg{ root }{ The root name of the files. 
Thus if nout=100 the files will be called names such as 
root001, root002 etc. If nout=1 then the output
is simply set = root.}
!!table

!!head2 Related commands

!!ref{tnadd.html}{tnadd}, !!ref{dvar.html}{dvar},
!!ref{dnadd.html}{dnadd}

!!end

 */

#include <climits>
#include <cstdlib>
#include <cmath>
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
    input.sign_in("seed",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nout",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("root",    Subs::Input::LOCAL,  Subs::Input::PROMPT);

    std::string intrail;
    input.get_value("trail", intrail, "trail", "input trailed spectrum");
    Subs::INT4 seed;
    input.get_value("seed",  seed, 554461, 1, INT_MAX, "random number seed");
    seed = -seed;
    int nout;
    input.get_value("nout",  nout, 10, 1, INT_MAX, "number of files to output");
    std::string root;
    input.get_value("root", root, "root", "root or base name for output files");

    Trail trail(intrail), tout;
    int tnum = 10, ndg=1;
    while(tnum <= nout){
      tnum *= 10;
      ndg++;
    }
    int nspec = trail.nspec(), npix=trail.npix();

    // array for counting choices

    int ntot = nspec*npix;
    int count[ntot], i, j, off;

    tout = trail;
    for(int n=0; n<nout; n++){

      // selection with replacement

      for(i=0; i<ntot; i++)
	count[i] = 0;

      for(i=0; i<ntot; i++){
	j = int(ntot*Subs::ran2(seed));
	if(j < 0 || j >= ntot)
	  throw Tomog::Tomog_Error("Bootstrap choice out of range!!");
	count[j]++;
      }
    
      // modify uncertainties
      for(int ns = 0; ns < nspec; ns++){
	off = npix*ns;
	for(int np = 0; np < npix; np++){
	  if(count[off+np]){
	    tout.error()[ns][np] = trail.error()[ns][np]/sqrt(float(count[off+np]));
	  }else{
	    tout.error()[ns][np] = -fabs(trail.error()[ns][np]);
	  }
	}
      }

      // write out file
      if(nout == 1){
	tout.write(root);
      }else{
	tout.write(root + Subs::str(int(n+1), ndg));
      }
    }
  }

  catch(const Trail::Trail_Error& err){
    std::cerr << "Trail::Trail_Error exception: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  catch(const Tomog::Tomog_Error& err){
    std::cerr << "Tomog::Tomog_Error exception: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}












