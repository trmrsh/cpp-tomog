/*

!!begin
!!title  Computes mean and rms scatter 
!!author T.R.Marsh
!!created  13 September 2000
!!revised  22 June 2003
!!root   dvar
!!index  dvar
!!descr  Computes mean and rms scatter of images
!!css   style.css
!!class  Noise
!!class  Doppler images
!!head1  dvar - computes mean and rms scatter of images

!!emph{dvar} calculates the mean and RMS scatter on every pixel of
a series. The results are output as images.

!!head2 Invocation

dvar list mean rms!!break

!!head2 Arguments

!!table
!!arg{ list  }{ list of images.}
!!arg{ mean  }{ mean of the images.}
!!arg{ rms   }{ rms of the images.}
!!table

!!end

*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <assert.h>
#include "trm_subs.h"
#include "trm_input.h"
#include "trm_tomog.h"
#include "trm_dmap.h"

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("list",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("mean",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("rms",    Subs::Input::GLOBAL, Subs::Input::PROMPT);

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

    std::string mean_name;
    input.get_value("mean",  mean_name, "mean", "output mean");
    std::string rms_name;
    input.get_value("rms",   rms_name, "rms", "output rms");

    Dmap mean(fname[0]), dummy;

    // First pass. Read in maps and compute their mean
    for(size_t nf=1; nf<fname.size(); nf++){
      dummy.read(fname[nf]);
      mean += dummy;
    }
    mean /= float(fname.size());
    mean.write(mean_name);

    // Second pass. Read in maps and compute their variance
    Dmap rms = mean;
    rms = 0;
    for(size_t nf=0; nf<fname.size(); nf++){
      dummy.read(fname[nf]);
      rms += Subs::sqr(dummy-mean);
    }
    rms /= float(std::max(size_t(1),fname.size()-1));
    rms.sqrt();
    rms.write(rms_name);
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




