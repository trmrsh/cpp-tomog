/*

!!begin
!!title  Contour plot a Doppler image
!!author T.R.Marsh
!!created  03 September 2000
!!revised  21 June 2003
!!root   dcont
!!index  dcont
!!descr  Contour plot a Doppler image
!!css   style.css
!!class  Display
!!class  Doppler images
!!head1  dcont - make a contour plot of a Doppler image.

dcont plots a Doppler image in contours of arbitrary level.

!!head2 Invocation

dcont input [device] nwave ngamma x1 x2 y1 y2 [title width csize] n*(lstyle clev)

!!head2 Arguments

!!table
!!arg{ input     }{ name of Doppler image file to plot.}
!!arg{ device    }{ the plot device.}
!!arg{ nwave     }{ the wavelength to plot.}
!!arg{ ngamma    }{ the systemic velocity to plot.}
!!arg{ x1        }{ lower x limit.}
!!arg{ x2        }{ upper x limit.}
!!arg{ y1        }{ lower y limit.}
!!arg{ y2        }{ upper y limit.}
!!arg{ title     }{ plot title.}
!!arg{ width     }{ width of plot; default=0.}
!!arg{ csize     }{ characer size of labels.}
!!arg{ lstyle    }{ line style for following contour level.}
!!arg{ contour   }{ contour level.}
!!table

!!end

*/

#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <vector>
#include "cpgplot.h"
#include "trm_subs.h"
#include "trm_constants.h"
#include "trm_plot.h"
#include "trm_input.h"
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
    input.sign_in("ngamma",  Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("x1",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("x2",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("y1",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("y2",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("title",   Subs::Input::LOCAL,  Subs::Input::NOPROMPT);
    input.sign_in("width",   Subs::Input::LOCAL,  Subs::Input::NOPROMPT);
    input.sign_in("csize",   Subs::Input::LOCAL,  Subs::Input::NOPROMPT);
    input.sign_in("lstyle",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("contour", Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("more",    Subs::Input::LOCAL,  Subs::Input::PROMPT);

    std::string infile;
    input.get_value("input",  infile, "input", "file to plot");
    std::string device;
    input.get_value("device",  device, "/xs", "plot device");
    Dmap map(infile);
    int nwave;
    input.get_value("nwave", nwave,   1, 1, map.nwave(),  "which wavelength to operate on");
    nwave--;
    int ngamma;
    input.get_value("ngamma", ngamma, 1, 1, map.ngamma(), "which systemic velocity to operate on");
    ngamma--;
    float vpix = map.vpix();
    int nx = map[nwave][ngamma].ncol(); 
    int ny = map[nwave][ngamma].nrow();
    
    float x1;
    input.get_value("x1", x1, -float(vpix*nx/2.), -FLT_MAX, FLT_MAX, "lower X limit (km/s)");
    float x2;
    input.get_value("x2", x2,  float(vpix*nx/2.), -FLT_MAX, FLT_MAX, "upper X limit (km/s)");
    float y1;
    input.get_value("y1", y1, -float(vpix*ny/2.), -FLT_MAX, FLT_MAX, "lower Y limit (km/s)");
    float y2;
    input.get_value("y2", y2,  float(vpix*ny/2.), -FLT_MAX, FLT_MAX, "upper Y limit (km/s)");
    std::string title;
    input.get_value("title",  title, "Doppler map", "plot title");
    float width;
    input.get_value("width", width, 0.f, 0.f, 100.f, "plot width (inches)");
    float csize;
    input.get_value("csize", csize, 0.f, 1.f, 100.f, "character size");

    // Read in line style and contour levels
    std::vector<int>   vstyle;
    std::vector<float> vcont;
    bool more = true;
    int style;
    float cont;
    while(more){
      input.get_value("lstyle",  style, 1, 1, 4, "line style"); 
      input.get_value("contour", cont,  0.f, -FLT_MAX, FLT_MAX, "contour level"); 
      vstyle.push_back(style);
      vcont.push_back(cont);
      input.get_value("more", more, true, "enter another style/contour pair?");
    } 

    // display image
    Subs::Plot plot(device);
    cpgpap(width,1.);
    cpgsch(csize);
    cpgscf(2);
    cpgvstd();
    cpgwnad(x1,x2,y1,y2);
    cpgsci(1);

    float tr[6];
    tr[0] = vpix*(-1.-(nx-1)/2.);
    tr[1] = vpix;
    tr[2] = 0.; 
    tr[3] = vpix*(-1.-(ny-1)/2.);
    tr[4] = 0.;
    tr[5] = vpix;

    float c[1];
    // Plot
    for(size_t i=0; i<vcont.size(); i++){
      cpgsls(vstyle[i]);
      c[0] = vcont[i];
      pgcont(map[nwave][ngamma],c,-1,tr);
    }

    cpgsci(4);
    cpgsls(1);
    cpgbox("BCNST",0.,0,"BCNST",0.,0);
    cpgsci(2);
    cpglab("V\\dx\\u (km s\\u-1\\d)","V\\dy\\u (km s\\u-1\\d)",title.c_str());
  }

  catch(const Dmap::Dmap_Error& err){
    std::cerr << "Dmap_Error: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  catch(const std::string& err){
    std::cerr << "Unrecognised std::string exception: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}
