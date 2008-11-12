/*

!!begin
!!title  Extracts symmetric part
!!author T.R.Marsh
!!created 13 September 2000
!!revised 22 June 2003
!!root   dsymm
!!index  dsymm
!!descr  Extracts azimuthally symmetric part of an image
!!css   style.css
!!class  Doppler images
!!head1  dsymm - extracts azimuthally symmetric part of an image

!!emph{dsymm} computes and returns the azimuthally symmetric part of an
image. The user must specify a centre of symmetry. The routine
works by interpolating the images onto circles centred on the
specified centre, sampling these, and then calculating the median
value. This median is then used to set the pixel values in the output.

The sampling is done in radial and azimuthal increments of 0.5 pixels.

!!head2 Invocation

dsymm map x0 y0 output!!break

!!head2 Arguments

!!table
!!arg{ map    }{ map to take symmetric part of.}
!!arg{ x0     }{ X position of centre of symmetry (km/s).}
!!arg{ y0     }{ Y position of centre of symmetry (km/s).}
!!arg{ output }{ output file name '-' for standard output.}
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
    input.sign_in("map",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("x0",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("y0",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("output",  Subs::Input::LOCAL, Subs::Input::PROMPT);

    std::string infile;
    input.get_value("map",  infile, "input", "file to take symmetric part of");
    float x0;
    input.get_value("x0", x0, 0.f, -FLT_MAX, FLT_MAX, "X centre of gaussian spot (km/s)");
    float y0;
    input.get_value("y0", y0, 0.f, -FLT_MAX, FLT_MAX, "Y centre of gaussian spot (km/s)");
    std::string outfile;
    input.get_value("output", outfile, "output", "output file");

    Dmap map(infile);

    const int MXBUFF = 5000;
    float r[MXBUFF], f[MXBUFF], t[MXBUFF];
    float vpix = map.vpix();

    int nx, ny, nrad, naz, nt1, ntot, na, ix, iy, j, nr;
    float dx, dy, rmax, x, y, xc, yc, fx, fy, rt, az;

    for(int nw=0; nw<map.nwave(); nw++){
      for(int ng=0; ng<map.ngamma(); ng++){
	
	// first compute radial profile
	ny   = map[nw][ng].nrow();
	nx   = map[nw][ng].ncol();
	xc   = x0/vpix + (nx-1)/2.; 
	yc   = y0/vpix + (ny-1)/2.; 
	
	dx   = std::max( fabs(vpix*(nx-1)/2.-x0), fabs(-vpix*(nx-1)/2.-x0));
	dy   = std::max( fabs(vpix*(ny-1)/2.-y0), fabs(-vpix*(ny-1)/2.-y0));
	rmax = sqrt(dx*dx+dy*dy)/vpix + 1.;
	nrad = int(2.*rmax);

	for(nr=0, ntot=0; nr<nrad; nr++){
	  r[nr] = rmax*nr/(nrad-1);
	  naz   = std::max(1,int(2.*Constants::TWOPI*r[nr]));
	  for(na=0, nt1=0; na<naz; na++){
	    az = Constants::TWOPI*na/naz;
	    x  = xc + r[nr]*cos(az);
	    y  = yc + r[nr]*sin(az);
	    if(x > 0.f && x < float(nx-1) && y > 0.f && y < float(ny-1)){
	      ix = int(x); iy = int(y);
	      fx = x - ix; fy = y - iy;
	      t[nt1++] = 
		(1.f-fx)*(1.f-fy)*map[nw][ng][iy][ix] + fx*(1.f-fy)*map[nw][ng][iy][ix+1] +
		(1.f-fx)*fy*map[nw][ng][iy+1][ix]     + fx*fy*map[nw][ng][iy+1][ix+1];
	    }
	  }
	  ntot = nr;
	  if(nt1 == 0) break;

	  Subs::quicksort(t,nt1);
	  f[nr] = t[nt1/2];
	  if(nt1 % 2 == 0) f[nr] = (f[nr] + t[nt1/2+1])/2.;
	}

	// now translate to output
	for(iy=0; iy<ny; iy++){
	  dy  = Subs::sqr(iy - yc);
	  for(ix=0; ix<nx; ix++){
	    rt  = sqrt(Subs::sqr(ix-xc)+dy);
	    j   = Subs::locate(r,ntot,rt);

	    if(j > 0 && j < ntot){
	      map[nw][ng][iy][ix] = ((r[j]-rt)*f[j-1]+(rt-r[j-1])*f[j])/(r[j]-r[j-1]);
	    }else{
	      map[nw][ng][iy][ix] = 0.f;
	    }
	  }
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

  exit(EXIT_SUCCESS);
}



