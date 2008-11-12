/*

!!begin
!!title  Plots distribution of flux in gamma
!!author T.R.Marsh
!!date   10 April 2001
!!root   dgdist
!!index  dgdist
!!descr  Plots distribution of flux in gamma
!!css   style.css
!!class  Doppler images
!!head1  dgdist - plots distribution of flux in gamma

!!emph{dgdist} plots distributin of flux in gamma. 

!!head2 Invocation

dgdist input [device] x1 x2 y1 y2 [title width csize]!!break

!!head2 Arguments

!!table
!!arg{ input  }{ input image ('-' for standard input.}
!!arg{ device }{ plot device}
!!arg{ nwave  }{ wavelength }
!!arg{ x1     }{ left-hand x limit of gamma}
!!arg{ x2     }{ right-hand x limit of gamma}
!!arg{ y1     }{ lower limit of plot}
!!arg{ y2     }{ upper limit of plot}
!!arg{ title  }{ title of plot}
!!arg{ width  }{ plot width, inches (default=device default)}
!!arg{ csize  }{ character size (default=1.5)}
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
#include "trm_plot.h"
#include "trm_tomog.h"
#include "trm_dmap.h"

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("input",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("device",  Subs::Input::GLOBAL, Subs::Input::NOPROMPT);
    input.sign_in("nwave",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("x1",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("x2",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("y1",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("y2",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("title",   Subs::Input::LOCAL,  Subs::Input::NOPROMPT);
    input.sign_in("width",   Subs::Input::LOCAL,  Subs::Input::NOPROMPT);
    input.sign_in("csize",   Subs::Input::LOCAL,  Subs::Input::NOPROMPT);

    std::string infile;
    input.get_value("input",  infile, "input", "file to plot");
    std::string device;
    input.get_value("device",  device, "/xs", "plot device");
    Dmap map(infile);
    int nwave;
    input.get_value("nwave", nwave, 1, 1, map.nwave(),  "which wavelength to operate on");
    nwave--;
    
    float x1;
    input.get_value("x1", x1, map.gamma(0), -FLT_MAX, FLT_MAX, "lower X limit (km/s)");
    float x2;
    input.get_value("x2", x2, map.gamma(map.ngamma()-1), -FLT_MAX, FLT_MAX, "upper X limit (km/s)");

    Subs::Array1D<float> flux(map.ngamma());
    for(int i=0; i<map.ngamma(); i++)
      flux[i] = sum(map[nwave][i]);

    float y1;
    input.get_value("y1", y1, 0.f, -FLT_MAX, FLT_MAX, "lower Y limit (km/s)");
    float y2;
    input.get_value("y2", y2, 1.1f*max(flux), -FLT_MAX, FLT_MAX, "upper Y limit (km/s)");
    std::string title;
    input.get_value("title",  title, "Doppler map", "plot title");
    float width;
    input.get_value("width", width, 0.f, 0.f, 100.f, "plot width (inches)");
    float csize;
    input.get_value("csize", csize, 0.f, 1.f, 100.f, "character size");

    // Plot
    Subs::Plot plot(device);
    cpgpap(width,1.);
    cpgsch(csize);
    cpgscf(2);
    cpgsci(4);
    cpgenv(x1,x2,y1,y2,0,0);
    cpgsci(2);
    cpglab("Systemic velocity, \\gg (km s\\u-1\\d)","Flux",title.c_str());
    cpgsci(1);
    pgline(map.gamma(),flux);
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
    std::cerr << "Memory allocation error" << std::endl;
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}


