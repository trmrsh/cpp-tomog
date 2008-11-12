/*

!!begin
!!title  Compute default image
!!author T.R.Marsh
!!date   5 April 2001
!!root   ddef
!!index  ddef
!!descr  Computes default image
!!css   style.css
!!class  Doppler images
!!head1  ddef - computes gaussian default image

!!emph{ddef} adds a gaussian spot to an image. 

!!head2 Invocation

dspot input nwave x0 y0 height width gamma wgamma output!!break

!!head2 Arguments

!!table
!!arg{ input  }{ input image (e.g. from !!emph{dinit}, '-' for standard input.}
!!arg{ fwhm   }{ image fwhm, pixels}
!!arg{ gfwhm  }{ gamma fwhm, pixels}
!!arg{ output }{ output file name '-' for standard output.}
!!table

!!end

*/


#include <stdlib.h>
#include <iostream.h>
#include <math.h>
#include <assert.h>
#include "trm_subs.h"
#include "trm_dmap.h"
#include "trm_constants.h"
#include "trm_buffer2d.h"

int main(int argc, char* argv[]){

  // -h for help, anything else for one liner

  if(argc == 2 && strcmp(argv[1],"-h") == 0){
    std::cerr << "\nddef computes default image\n\n";
    std::cerr << "usage: ddef input fwhm gfwhm output\n\n";
    std::cerr <<  "  input  -- input image\n";
    std::cerr <<  "  fwhm   -- image fwhm pixels.\n";
    std::cerr <<  "  gfwhm  -- gamma fwhm pixels.\n";
    std::cerr <<  "  output -- output file name ('-' for standard output)\n\n";
    exit(EXIT_FAILURE);
  }else if(argc != 10){
    std::cerr << "usage: dspot input nwave x0 y0 height width gamma wgamma output\n\n";
    exit(EXIT_FAILURE);
  }

  // automatic counting, just need to get order right.
  // a few assertions are checked.
 
  try{

    int narg = 0;

    Dmap map(argv[++narg]);
    size_t nwave = atoi(argv[++narg]); assert(nwave > 0 && nwave <= map.nwave());
    size_t nside = map.nside();
    float vpix   = map.get_vpix();
    nwave--;

    float x       = atof(argv[++narg]);
    float y       = atof(argv[++narg]);
    float h       = atof(argv[++narg]);
    float w       = atof(argv[++narg]); assert(w > 0.);
    float gamma   = atof(argv[++narg]);
    float wgamma  = atof(argv[++narg]); assert(wgamma > 0.);

    std::string output(argv[++narg]);

    Buffer2D<float> spot(nside,nside);
    float efac = 1./sqr(w/vpix/EFAC)/2.;
    float cen  = float(nside-1)/2.;
    x = x/vpix + cen;
    y = y/vpix + cen;
    for(size_t iy=0; iy<nside; iy++)
      for(size_t ix=0; ix<nside; ix++)
	spot[iy][ix] = h*exp(-efac*(sqr(ix-x)+sqr(iy-y)));

    for(size_t ig=0; ig<map.ngamma(); ig++)
      map[nwave][ig] += float(exp(-sqr((map.get_gamma(ig)-gamma)/
				       (wgamma/EFAC))/2.))*spot;
    map.write(output);
  }

  catch(Dmap::Dmap_error e){
    std::cerr << "Dmap:: error!!" << std::endl;
    std::cerr << e.get_mess() << std::endl;
    exit(EXIT_FAILURE);
  }
  catch(bad_alloc){
    std::cerr << "Memory allocation error" << std::endl;
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}


