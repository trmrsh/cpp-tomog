/*

!!begin
!!title  Rank pixel-by-pixel
!!author T.R.Marsh
!!date   13 September 2000
!!date   22 June 2003
!!root   drank
!!index  drank
!!descr  Rank an image pixel-by-pixel relative to others
!!css   style.css
!!class  Arithematic
!!class  Doppler images
!!head1  drank - rank an image pixel-by-pixel relative to others

!!emph{drank} takes one image and then works out where each of its pixels
ranks in comparison to a set of other images contained in a list.
It is useful in monte carlo runs of significance.

Each pixel of the output image contains the number of images that the 
pixel in question was larger than. Thus if the list contains 10 images,
a value of 0 would indicate that the pixel was less than its equivalent
in the other 10, while a value of 10 would show that it was higher
than all of them.

!!head2 Invocation

drank input list output!!break

!!head2 Arguments

!!table
!!arg{ input  }{ name of Doppler image file.}
!!arg{ list   }{ list of Doppler images for reference.}
!!arg{ output }{ output Doppler image}
!!table

!!end

*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include "trm_subs.h"
#include "trm_input.h"
#include "trm_tomog.h"
#include "trm_dmap.h"

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("input",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("list",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("output",  Subs::Input::GLOBAL, Subs::Input::PROMPT);

    std::string infile;
    input.get_value("input",  infile, "input", "file to plot");

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

    std::string outfile;
    input.get_value("output", outfile, "output", "output file");

    Dmap dmap(infile);
    Dmap dummy, output = dmap;
    output = 0.;

    // First pass. Read in maps and compute their mean
    for(size_t nf=0; nf<fname.size(); nf++){
      dummy.read(fname[nf]);
      if(!match(dummy,dmap))
	throw Tomog::Input_Error("Non-matching image.");

      for(int nw=0; nw<dmap.nwave(); nw++){
	for(int ng=0; ng<dmap.ngamma(); ng++){
	  for(int iy=0; iy<dmap[nw][ng].nrow(); iy++){
	    for(int ix=0; ix<dmap[nw][ng].ncol(); ix++){
	      if(dmap[nw][ng][iy][ix] > dummy[nw][ng][iy][ix]) output[nw][ng][iy][ix]++;
	    }
	  }
	}
      }
    }
    output /= fname.size();
    output.write(outfile);
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


