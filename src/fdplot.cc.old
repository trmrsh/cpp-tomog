/* 

!!begin
!!title  Fancy plot of Doppler images as greyscales
!!author T.R.Marsh
!!created 13 September 2000
!!revised 22 June 2003
!!root   fdplot
!!index  fdplot
!!descr  Fancy plot of Doppler images as greyscale images
!!css   style.css
!!class  Display
!!class  Doppler images
!!head1  fdplot - fancy plot of Doppler images as greyscale images

!!emph{fdplot} is a fancy version of !!ref{dplot.html}{dplot}
and allows one to plot multiple Doppler images butted
together.

!!head2 Invocation

fdplot [device] nx ny x1 x2 y1 y2 width csize reverse (map nwave ngamma low high title more)*n!!break

Arguments in brackets can be repeated many times.


!!head2 Arguments

!!table
!!arg{ device }{ plot device}
!!arg{ nx     }{ number of panels in x}
!!arg{ ny     }{ number of panels in y}
!!arg{ x1     }{ lower x limit}
!!arg{ x2     }{ upper x limit}
!!arg{ y1     }{ lower y limit}
!!arg{ y2     }{ upper y limit}
!!arg{ width  }{ plot width, inches}
!!arg{ csize  }{ character size}
!!arg{ reverse}{ reverse usual background colour (i.e. white becomes black and vice versa)}
!!arg{ map    }{ image to plot}
!!arg{ nwave  }{ wavelength number of image}
!!arg{ ngamma }{ systemic velocity number of image to plot}
!!arg{ low    }{ sets lower plot level}
!!arg{ high   }{ sets upper plot level}
!!arg{ title  }{ plot title}
!!arg{ more   }{ controls whether you want to input another file}
!!table

The last set can be done for every image.

!!end

*/

#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <string>
#include "cpgplot.h"
#include "trm_subs.h"
#include "trm_plot.h"
#include "trm_input.h"
#include "trm_tomog.h"
#include "trm_dmap.h"

// Structure for each panel
struct Panel {

  Panel(int nside, std::string title_, float low_, float high_, float vpix_) 
    : dmap(nside,nside), title(title_), low(low_), high(high_), vpix(vpix_) {}

  Subs::Array2D<float> dmap;
  std::string title;
  float low, high, vpix;
};

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("device",  Subs::Input::GLOBAL, Subs::Input::NOPROMPT);
    input.sign_in("nx",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("ny",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("x1",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("x2",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("y1",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("y2",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("width",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("csize",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("reverse", Subs::Input::LOCAL,  Subs::Input::NOPROMPT);
    input.sign_in("map",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nwave",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("ngamma",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("low",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("high",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("title",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("more",    Subs::Input::LOCAL,  Subs::Input::PROMPT);

    std::string device;
    input.get_value("device",  device, "/xs", "plot device");
    int nx;
    input.get_value("nx", nx, 3, 1, 20,  "number of panels in X");
    int ny;
    input.get_value("ny", ny, 3, 1, 20,  "number of panels in Y");
    float x1;
    input.get_value("x1", x1, -1000.f, -FLT_MAX, FLT_MAX, "lower X limit (km/s)");
    float x2;
    input.get_value("x2", x2,  1000.f, -FLT_MAX, FLT_MAX, "upper X limit (km/s)");
    float y1;
    input.get_value("y1", y1, -1000.f, -FLT_MAX, FLT_MAX, "lower Y limit (km/s)");
    float y2;
    input.get_value("y2", y2,  1000.f, -FLT_MAX, FLT_MAX, "upper Y limit (km/s)");
    float width;
    input.get_value("width", width, 0.f, 0.f, 100.f, "plot width (inches)");
    bool rev;
    input.get_value("reverse", rev, false, "reverse usual background colour?");
    float csize;
    input.get_value("csize", csize, 0.f, 1.f, 100.f, "character size");

    bool more = true;
    std::string infile, title;
    Dmap dmap;
    int nwave, ngamma;
    float low, high;
    std::vector<Panel> panel;
    Panel *pptr;
    while(more && int(panel.size()) < nx*ny){

      input.get_value("map",  infile, "input", "file to plot");
      dmap = Dmap(infile);
      input.get_value("nwave",  nwave, 1, 1, dmap.nwave(), "which wavelength to operate on");
      input.get_value("ngamma", ngamma, 1, 1, dmap.ngamma(), "which systemic velocity to operate on");
      nwave--; ngamma--;
      input.get_value("low",    low,  0.f,   -FLT_MAX, FLT_MAX, "lower intensity limit");
      input.get_value("high",   high, 100.f, -FLT_MAX, FLT_MAX, "upper intensity limit");
      input.get_value("title",  title, "title", "plot title");

      panel.push_back(Panel(dmap.nside(), title, low, high, dmap.vpix()));

      pptr = &panel[panel.size()-1];

      // transfer data to avoid trying to store too much data
      for(int iy=0; iy<dmap.nside(); iy++){
	for(int ix=0; ix<dmap.nside(); ix++){
	  pptr->dmap[iy][ix] = dmap[nwave][ngamma][iy][ix];
	}
      }

      if(int(panel.size()) < nx*ny) input.get_value("more", more, true, "next panel?");
    }

    // display images
    float le = 4., re = 3., be = 4., te = 4.;

    Subs::Plot plot(device);
    float xh, yh, aspect;
    float xr = nx*(x2-x1), yr = ny*(y2-y1);
    xh = xr/(1.-(le+re)*csize/40.);
    yh = yr/(1.-(be+te)*csize/40.);
    if(xh < yh){
      yh = yr + (be+te)*csize*xh/40.;
    }else{
      xh = xr + (le+re)*csize*yh/40.;
    }
    aspect = yh/xh;

    if(rev){
      cpgscr(0, 1., 1., 1.);
      cpgscr(1, 0., 0., 0.);
    }
    cpgpap(width,aspect);
    cpgsch(csize);
    cpgscf(2);
    float cs  = csize*std::min(1.f,aspect)/40.;
    float xv1 = le*cs;
    float xv2 = 1.-re*cs;
    float yv1 = be*cs/aspect;
    float yv2 = 1.-te*cs/aspect;

    float xtv1, xtv2, ytv1, ytv2;
    float vpix, tr[6];
    std::string xlab, ylab;
    char xopt[6], yopt[6];

    for(int iy=ny-1, npanel=0; iy>=0; iy--){
      ytv1 = yv1 + iy*(yv2-yv1)/ny;
      ytv2 = yv1 + (iy+1)*(yv2-yv1)/ny;
      for(int ix=0; ix<nx; ix++, npanel++){
	xtv1 = xv1 + ix*(xv2-xv1)/nx;
	xtv2 = xv1 + (ix+1)*(xv2-xv1)/nx;
	cpgsvp(xtv1,xtv2,ytv1,ytv2);
	cpgswin(x1,x2,y1,y2);
	cpgsci(1);
	
	if(npanel < int(panel.size())){

	  pptr = &panel[npanel];
	  int nside = pptr->dmap.ncol();
	  vpix  = pptr->vpix;
	  tr[0] = vpix*(-1.-(nside-1)/2.);
	  tr[1] = vpix;
	  tr[2] = 0.; 
	  tr[3] = vpix*(-1.-(nside-1)/2.);
	  tr[4] = 0.;
	  tr[5] = vpix;
	  pggray(pptr->dmap,pptr->high,pptr->low,tr);
	  cpgsci(4);
	  if(ix == 0){
	    strcpy(yopt,"BCNST");
	  }else{
	    strcpy(yopt,"BCST");
	  }
	  if(iy == 0){
	    strcpy(xopt,"BCNST");
	  }else{
	    strcpy(xopt,"BCST");
	  }
	  cpgbox(xopt,0.,0,yopt,0.,0);
	  cpgsci(2);
	  if(ix == 0){
	    ylab = "V\\dy\\u (km s\\u-1\\d)";
	  }else{
	    ylab = "";
	  }
	  if(iy == 0){
	    xlab = "V\\dx\\u (km s\\u-1\\d)";
	  }else{
	    xlab = "";
	  }
	  if(iy == ny-1){
	    cpglab(xlab.c_str(), ylab.c_str(), title.c_str());
	  }else{
	    cpglab(xlab.c_str(),ylab.c_str()," ");
	  }
	}
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
