/*

!!begin
!!title   Returns information on Doppler images and trailed spectra
!!author  T.R.Marsh
!!created 13 September 2000
!!revised 22 June 2003
!!root    dtinfo
!!index   dtinfo
!!descr   Returns informtion on Doppler images and trailed spectra
!!css   style.css
!!class   Doppler images
!!class   Trailed spectra
!!head1   dtinfo - returns information on Doppler images and trailed spectra

!!emph{dtinfo} comes back with basic information on Doppler images and
trailed spectra such as the number of images, their central wavelengths 
etc. It can recognise the difference between the two types and act accordingly.

!!head2 Invocation

dtinfo input!!break

!!head2 Arguments

!!table
!!arg{ input }{ name of Doppler image or trailed spectrum.}
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
#include "trm_trail.h"

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("input",  Subs::Input::GLOBAL, Subs::Input::PROMPT);

    std::string infile;
    input.get_value("input",  infile, "input", "file to add spot to");
    std::ifstream file(infile.c_str());
    if(!file) throw Tomog::Input_Error("Could not open file = " + infile);
    int tflag;
    file.read((char*)&tflag,sizeof(tflag));
    file.close();

    if(tflag == Dmap::flag){
      std::cout << "\n         File type: Doppler image\n";

      Dmap map(infile);
      std::cout << " Number of wavelengths = " << map.nwave() << std::endl;
      std::cout << " Number of gamma vels  = " << map.ngamma() << std::endl;
      std::cout << "Systemic velocity range = " << map.gamma(0) << " to " << map.gamma(map.ngamma()-1) << " km/s" << std::endl;
      std::cout << "   Velocity/pixel = " << map.vpix() << " km/s" << std::endl;
      std::cout << "      Pixels/side = " << map.nside() << std::endl; 
      std::cout << "Range: " << min(map) << " to " << max(map) << std::endl;
      for(int n=0; n<map.nwave(); n++)
	std::cout << "Wavelength number " << n+1 << " = " << map.wzero(n) << std::endl;

    }else if(tflag == Trail::flag){
      std::cout << "\n                 File type: Trailed spectrum\n";

      Trail trail(infile);
      std::cout << "Number of pixels/spectrum = " << trail.npix()      << std::endl;
      std::cout << "        Number of spectra = " << trail.nspec()     << std::endl;
      std::cout << "           Velocity/pixel = " << trail.vpix()      << " km/s" << std::endl;
      std::cout << "       Central wavelength = " << trail.wzero()     << std::endl;
      std::cout << "                     Range: " << min(trail.data()) << " to " << max(trail.data()) << std::endl;
    }else{
      std::cout << "\nFile type not recognised!\n";
    }
  }

  catch(const Dmap::Dmap_Error& err){
    std::cerr << "Dmap::Dmap_Error exception: " << err << std::endl;
    exit(EXIT_FAILURE);
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


