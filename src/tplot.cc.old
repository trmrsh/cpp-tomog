/* 

!!begin
!!title  Plotting Doppler images as greyscales
!!author T.R.Marsh
!!date   13 September 2000
!!date   23 June 2003
!!root   tplot
!!index  tplot
!!descr  Plots a trailed spectrum
!!css   style.css
!!class  Display
!!class  Trailed spectra
!!head1  tplot - plots trailed spectra as greyscale images

!!emph{tplot} is the workhorse plot program for trailed spectra
within the C++ !!ref{index.html}{Doppler tomography} suite of programs.

!!head2 Invocation

tplot file [device] x1 x2 y1 y2 low high [title width csize aspect]!!break

Arguments in brackets are optional. 

!!head2 Arguments

!!table
!!arg{ trail   }{ name of trailed spectrum file.}
!!arg{ device  }{ plotting device}
!!arg{ x1    }{ left-hand x limit }
!!arg{ x2    }{ right-hand x limit}
!!arg{ y1    }{ lower y limit }
!!arg{ y2    }{ upper y limit }
!!arg{ low    }{ lower plot level }
!!arg{ high    }{ upper plot level }
!!arg{ title }{ plot title }
!!arg{ width }{ plot width, inches }
!!arg{ csize }{ character size }
!!arg{ aspect }{ aspect ratio (height/width)}
!!table

!!end

*/

#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <string>
#include <iostream>
#include "cpgplot.h"
#include "trm_subs.h"
#include "trm_input.h"
#include "trm_plot.h"
#include "trm_tomog.h"
#include "trm_trail.h"

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("input",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("device",  Subs::Input::GLOBAL, Subs::Input::NOPROMPT);
    input.sign_in("x1",      Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("x2",      Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("y1",      Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("y2",      Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("low",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("high",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("title",   Subs::Input::LOCAL,  Subs::Input::NOPROMPT);
    input.sign_in("width",   Subs::Input::LOCAL,  Subs::Input::NOPROMPT);
    input.sign_in("csize",   Subs::Input::LOCAL,  Subs::Input::NOPROMPT);
    input.sign_in("aspect",  Subs::Input::LOCAL,  Subs::Input::NOPROMPT);

    std::string infile;
    input.get_value("input",  infile, "input", "file to plot");
    std::string device;
    input.get_value("device",  device, "/xs", "plot device");

    Trail trail(infile);
      
    int nspec  = trail.nspec(); 
    int npix   = trail.npix();
    float vpix = trail.vpix();
    
    float x1;
    input.get_value("x1", x1, float(vpix*(-0.5-(npix-1)/2.)), -FLT_MAX, FLT_MAX, "lower X limit (km/s)");
    float x2;
    input.get_value("x2", x2, float(vpix*(npix-0.5-(npix-1)/2.)), -FLT_MAX, FLT_MAX, "upper X limit (km/s)");
    float y1;
    input.get_value("y1", y1, 0.5f, -FLT_MAX, FLT_MAX, "lower Y limit (in terms of spectrum numbers)");
    float y2;
    input.get_value("y2", y2, float(nspec+0.5), -FLT_MAX, FLT_MAX, "upper Y limit (in terms of spectrum numbers)");
    float low;
    input.get_value("low",  low, 0.f, -FLT_MAX, FLT_MAX, "lower intensity limit");
    float high;
    input.get_value("high", high, 100.f, -FLT_MAX, FLT_MAX, "upper intensity limit");
    std::string title;
    input.get_value("title",  title, "Doppler map", "plot title");
    float width;
    input.get_value("width", width, 0.f, 0.f, 100.f, "plot width (inches)");
    float csize;
    input.get_value("csize", csize, 0.f, 1.f, 100.f, "character size");
    float aspect;
    input.get_value("aspect", aspect, 0.f, 0.f, 100.f, "apect ratio (height/width)");

    // display image
    Subs::Plot plot(device);

    float xs1, xs2, ys1, ys2;
    cpgqvsz(1,&xs1,&xs2,&ys1,&ys2);

    if(aspect == 0.) aspect = (ys2-ys1)/(xs2-xs1);
    cpgpap(width,aspect);

    cpgsch(csize);
    cpgscf(2);
    cpgvstd();
    cpgswin(x1,x2,y1,y2);
    cpgsci(1);
    float tr[6];
    tr[0] = vpix*(-1.-(npix-1)/2.);
    tr[1] = vpix;
    tr[2] = 0.; 
    tr[3] = 0.;
    tr[4] = 0.;
    tr[5] = 1.;
    pggray(trail.data(),high,low,tr);
    cpgsci(4);
    cpgbox("BCNST",0.,0,"BCNST",0.,0);
    cpgsci(2);
    cpglab("Velocity (km s\\u-1\\d)","Spectrum number",title.c_str());
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
