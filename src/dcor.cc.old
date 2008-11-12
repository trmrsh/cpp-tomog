/*

!!begin
!!title  Correlation
!!author T.R.Marsh
!!date   13 September 2000
!!root   dcor
!!index  dcor
!!descr  Computes correlation between one pixel and rest of image
!!css   style.css
!!class  Arithematic
!!class  Doppler images
!!head1  dcor - computes correlation between one pixel and rest of image

dcor takes a list of images and determines the correlation between
a particular pixel and the rest of the image. The result is itself
an image, peaking at a value of 1 on the selected pixel.


!!head2 Invocation

dcor list nx ny nim output!!break

!!head2 Arguments

!!table
!!arg{ list   }{ file containing a series of image names.}
!!arg{ nwave  }{ wavelength number to process}
!!arg{ ngamma }{ gamma number to process}
!!arg{ nx     }{ x coordinate of pixel (0 to 1 less than x dimension).}
!!arg{ ny     }{ y coordinate of pixel (0 to 1 less than y dimension).}
!!arg{ output }{ output Doppler image}
!!table

!!end

*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "trm_subs.h"
#include "trm_tomog.h"
#include "trm_input.h"
#include "trm_dmap.h"

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("list",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("nwave",  Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("ngamma", Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("nx",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("ny",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("output", Subs::Input::LOCAL,  Subs::Input::PROMPT);

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
      
    // Read in first image
    Dmap mean(fname[0]);

    int nwave;
    input.get_value("nwave", nwave,   1, 1, mean.nwave(),  "which wavelength to operate on");
    nwave--;
    int ngamma;
    input.get_value("ngamma", ngamma, 1, 1, mean.ngamma(), "which systemic velocity to operate on");
    ngamma--;

    int nx;
    input.get_value("nx", nx, 1, 1, mean.nside(), "X pixel to correlate against");
    nx--;
    int ny;
    input.get_value("ny", ny, 1, 1, mean.nside(), "Y pixel to correlate against");
    ny--;

    std::string outfile;
    input.get_value("output", outfile, "output", "output file");

    Dmap dummy;
    for(size_t nf=1; nf<fname.size(); nf++){
      dummy.read(fname[nf]);
      mean += dummy;
    }
    mean /= float(fname.size());

    float mpix = mean[nwave][ngamma][ny][nx];

    // Second pass. Read in maps and compute correlations
    Dmap  var, cvar;
    var  = 0;
    cvar = 0;
    for(size_t nf=0; nf<fname.size(); nf++){
      dummy.read(fname[nf]);
      var  += Subs::sqr(dummy-mean);
      cvar += (dummy[nwave][ngamma][ny][nx]-mpix)*(dummy-mean);
    }
    var.sqrt();
    cvar /= var;
    cvar /= var[nwave][ngamma][ny][nx];
    cvar.write(outfile);
  }

  catch(const Dmap::Dmap_Error& err){
    std::cerr << "Dmap_Error: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  catch(const std::string& err){
    std::cerr << "string error: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}


