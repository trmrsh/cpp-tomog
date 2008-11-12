/* 

!!begin
!!title  Plotting Doppler images as greyscales
!!author T.R.Marsh
!!date   04 July 2000
!!root   dplot
!!index  dplot
!!descr  Plots a Doppler image as a greyscale image
!!css   style.css
!!class  Display
!!class  Doppler images
!!head1  dplot - plots Doppler images as greyscale images

!!emph{dplot} is the workhorse plot program for Doppler images
within the C++ !!ref{index.html}{Doppler 
tomography} suite of programs.

!!head2 Invocation

dplot input [device] nwave ngamma x1 x2 y1 y2 low high [title width csize]!!break

!!head2 Arguments

!!table
!!arg{ map     }{ name of Doppler image file.}
!!arg{ device  }{ plotting device}
!!arg{ nwave   }{ wavelength number to plot.}
!!arg{ ngamma  }{ gamma velocity slice to plot, 0 for sum of them all.}
!!arg{ x1      }{ left-hand x limit}
!!arg{ x2      }{ right-hand x limit}
!!arg{ y1      }{ lower y limit}
!!arg{ y2      }{ upper y limit}
!!arg{ low     }{ lower plot level}
!!arg{ high    }{ upper plot level}
!!arg{ title   }{ plot title}
!!arg{ width   }{ plot width, inches}
!!arg{ csize   }{ character size}
!!table

!!end

*/

#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <string>
#include "cpgplot.h"
#include "trm_subs.h"
#include "trm_array2d.h"
#include "trm_plot.h"
#include "trm_input.h"
#include "trm_tomog.h"
#include "trm_dmap.h"

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("map",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("device",  Subs::Input::GLOBAL, Subs::Input::NOPROMPT);
    input.sign_in("nwave",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("ngamma",  Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("x1",      Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("x2",      Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("y1",      Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("y2",      Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("low",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("high",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("title",   Subs::Input::LOCAL,  Subs::Input::NOPROMPT);
    input.sign_in("width",   Subs::Input::LOCAL,  Subs::Input::NOPROMPT);
    input.sign_in("csize",   Subs::Input::LOCAL,  Subs::Input::NOPROMPT);

    std::string infile;
    input.get_value("map",   infile, "input", "file to plot");
    std::string device;
    input.get_value("device",  device, "/xs", "plot device");
    Dmap dmap(infile);
    int nwave;
    input.get_value("nwave", nwave,   1, 1, dmap.nwave(),  "which wavelength to operate on");
    nwave--;
    int ngamma;
    bool sum;
    if(dmap.ngamma() > 1){
      input.get_value("ngamma", ngamma, 1, 0, dmap.ngamma(), "which systemic velocity to operate on");
      if(ngamma){
	ngamma--;
	sum = false;
      }else{
	sum = true;
      }
    }else{
      sum = false;
      ngamma = 0;
    }
      
    int nside = dmap.nside();
    float vpix   = dmap.vpix();
    
    float x1;
    input.get_value("x1", x1, -float(vpix*nside/2.), -FLT_MAX, FLT_MAX, "lower X limit (km/s)");
    float x2;
    input.get_value("x2", x2,  float(vpix*nside/2.), -FLT_MAX, FLT_MAX, "upper X limit (km/s)");
    float y1;
    input.get_value("y1", y1, -float(vpix*nside/2.), -FLT_MAX, FLT_MAX, "lower Y limit (km/s)");
    float y2;
    input.get_value("y2", y2,  float(vpix*nside/2.), -FLT_MAX, FLT_MAX, "upper Y limit (km/s)");
    float low;
    input.get_value("low",  low, 0.f, -FLT_MAX, FLT_MAX, "lower intensity limit");
    float high;
    input.get_value("high", high, dmap.max(), -FLT_MAX, FLT_MAX, "upper intensity limit");

    std::string title;
    input.get_value("title",  title, "Doppler map", "plot title");
    float width;
    input.get_value("width", width, 0.f, 0.f, 100.f, "plot width (inches)");
    float csize;
    input.get_value("csize", csize, 0.f, 1.f, 100.f, "character size");

    // display image
    Subs::Plot plot(device);
    cpgpap(width,1.);
    cpgsch(csize);
    cpgscf(2);
    cpgvstd();
    cpgwnad(x1,x2,y1,y2);
    cpgsci(1);

    float tr[6];
    tr[0] = vpix*(-1.-(nside-1)/2.);
    tr[1] = vpix;
    tr[2] = 0.; 
    tr[3] = vpix*(-1.-(nside-1)/2.);
    tr[4] = 0.;
    tr[5] = vpix;

    if(sum){
      Subs::Array2D<float> temp = dmap[nwave][0];
      for(int i=1; i<dmap.ngamma(); i++)
	temp += dmap[nwave][i];
      pggray(temp,high,low,tr);
    }else{
      pggray(dmap[nwave][ngamma],high,low,tr);
    }
    cpgsci(4);
    cpgbox("BCNST",0.,0,"BCNST",0.,0);
    cpgsci(2);
    cpglab("V\\dx\\u (km s\\u-1\\d)","V\\dy\\u (km s\\u-1\\d)",title.c_str());
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
