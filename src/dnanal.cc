/*

!!begin
!!title  Noise analysis of Doppler images
!!author T.R.Marsh
!!date   13 September 2000
!!date   22 June 2003
!!root   dnanal
!!index  dnanal
!!descr  Noise analysis of Doppler images
!!css   style.css
!!class  Noise
!!class  Doppler images
!!head1  dnanal - carries out noise analysis of Doppler images

!!emph{dnanal} works out the RMS variation in circular apertures covering
a range of radii given a list of Doppler images. The circular
apertures are centred on a particular x,y position and a series
of radii is taken from a specified inner to outer value. The results
are sent to standard output.

!!head2 Invocation

dnanal list x y nim r1 r2 nr!!break

!!head2 Arguments

!!table
!!arg{ list   }{ list of Doppler images.}
!!arg{ x      }{ x centre of circles (pixels).}
!!arg{ y      }{ y centre of circles (pixels).}
!!arg{ nim    }{ image number to analyse.}
!!arg{ r1     }{ first radius of circles in pixels.}
!!arg{ r2     }{ last radius of circles in pixels.}
!!arg{ nrad   }{ number of radii.}
!!table

!!end

*/

#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <fstream>
#include "trm_subs.h"
#include "trm_buffer2d.h"
#include "trm_array2d.h"
#include "trm_input.h"
#include "trm_tomog.h"
#include "trm_dmap.h"

void tcirc(const Subs::Buffer2D<float>& img, float x, float y, float r, float& sum, int& npix);

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("list",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("x",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("y",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nwave",  Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("ngamma", Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("r1",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("r2",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nrad",   Subs::Input::LOCAL,  Subs::Input::PROMPT);

    // Read in list of file names
    std::string flist;
    input.get_value("list",  flist, "flist", "list of maps");
    std::ifstream list(flist.c_str());
    if(!list)
      throw Tomog::Input_Error("Could not open " + flist);
    std::string file;
    std::vector<std::string> fname;
    while(list >> file){
      fname.push_back(file);
    }
    list.close();
    if(fname.size() == 0)
      throw Tomog::Input_Error("No file names loaded from " + flist);

    float x;
    input.get_value("x", x, 0.f, -FLT_MAX, FLT_MAX, "X centre of computations (pixels)");
    float y;
    input.get_value("y", y, 0.f, -FLT_MAX, FLT_MAX, "Y centre of computations (pixels)");
    x--; y--;

    Dmap mean(fname[0]);

    int nwave;
    input.get_value("nwave", nwave,   1, 1, mean.nwave(),  "which wavelength to operate on");
    nwave--;
    int ngamma;
    input.get_value("ngamma", ngamma, 1, 1, mean.ngamma(), "which systemic velocity to operate on");
    ngamma--;

    float r1;
    input.get_value("r1", r1, 0.f, 0.f, FLT_MAX, "inner radius (pixels)");
    float r2;
    input.get_value("r2", r2, std::max(r1, 100.f), r1, FLT_MAX, "outer radius (pixels)");
    int nrad;
    input.get_value("nrad", nrad, 10, 2, 10000, "number of radii");

    Dmap dummy;
    for(size_t nf=1; nf<fname.size(); nf++){
      dummy.read(fname[nf]);
      mean += dummy;
    }
    mean /= float(fname.size());

    // compute means over circles

    std::vector<float> scirc(nrad), rms(nrad);
    std::vector<int> ncirc(nrad);
    float r;
    for(int n=0; n<nrad; n++){
      r = r1 + (r2-r1)*n/(nrad-1);
      tcirc(mean[nwave][ngamma], x, y, r, scirc[n], ncirc[n]);
    }

    // Second pass. Read in maps and compute rms
    float sum;
    int npix;
    for(size_t nf=0; nf<fname.size(); nf++){
      dummy.read(fname[nf]);
      for(int n=0; n<nrad; n++){
	if(ncirc[n]){	
	  r = r1 + (r2-r1)*n/(nrad-1);
	  tcirc(dummy[nwave][ngamma], x, y, r, sum, npix);
	  rms[n] += Subs::sqr(sum-scirc[n]);
	}
      }
    }
    
    // write out results
    for(int n=0; n<nrad; n++){
      if(ncirc[n]){	
	r = r1 + (r2-r1)*n/(nrad-1);
	std::cout << r << " " << ncirc[n] << " " << scirc[n] << " " << sqrt(rms[n]/std::max(int(fname.size()-1),1)) << std::endl;
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

// sums up flux over circle centred x,y, radius r

void tcirc(const Subs::Buffer2D<float>& img, float x, float y, float r, float& sum, int& npix){
  int ny = img.nrow();
  int nx = img.ncol();
  float rsq = r*r;

  sum  = 0.;
  npix = 0;
  float extra;
  for(int iy=0; iy<ny; iy++){
    extra = Subs::sqr(iy-y);
    if(extra < rsq){
      for(int ix=0; ix<nx; ix++){
	if((extra + Subs::sqr(ix-x)) < rsq){
	  sum += img[iy][ix];
	  npix++;
	}
      }
    }
  }
}

