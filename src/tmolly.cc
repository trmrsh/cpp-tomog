/*

!!begin
!!title  Converts a trail to molly format
!!author T.R.Marsh
!!date   30 September 2003
!!root   tmolly
!!index  tmolly
!!descr  Converts a trail to molly format
!!css   style.css
!!class  Trailed spectra
!!head1  tmolly - converts a trail to molly format

!!emph{tmolly} reads in atrail and outputs an equivalent molly file.

!!head2 Invocation

tmolly trail molly!!break

!!head2 Arguments

!!table
!!arg{trail }{ name of input trailed spectrum.}
!!arg{molly }{ name of output molly file.}
!!table

!!end

 */

#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <string>
#include <iostream>
#include "trm_subs.h"
#include "trm_constants.h"
#include "trm_input.h"
#include "trm_colly.h"
#include "trm_tomog.h"
#include "trm_trail.h"

void strint(const unsigned int num, const unsigned int nd, char* intstr);

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("trail",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("molly",   Subs::Input::LOCAL,  Subs::Input::PROMPT);

    std::string infile;
    input.get_value("trail",  infile, "input", "input trailed spectrum");
    std::string outfile;
    input.get_value("molly",  outfile, "input", "input trailed spectrum");

    Trail trail(infile);

    // Generate suitable header 
    Subs::Header head;
    const int fcode = 5;
    const std::string units = "MILLIJANSKYS    ";
    head.set("Xtra",       new Subs::Hdirectory("Extra information on spectrum"));
    head.set("Xtra.FCODE", new Subs::Hint(fcode,"molly format code"));
    head.set("Xtra.UNITS", new Subs::Hstring(units, "Flux units"));
    head.set("Xtra.NPIX",  new Subs::Hint(trail.npix(),"Number of pixels"));
    head.set("Xtra.NARC",  new Subs::Hint(-2, "Number of arc coefficients"));    
    std::vector<double> arc(2);
    arc[1] = trail.npix()*trail.vpix()/(Constants::C/1000.);
    arc[0] = log(trail.wzero()) - arc[1]*(trail.npix()+1)/(2.*trail.npix());
    head.set("Xtra.ARC",  new Subs::Hdvector(arc, "Arc coefficients"));    
    head.set("Object",    new Subs::Hstring("Doppler trail spectrum", "Object name"));

    // Previous parts are the same for each spectrum. Now open output file and
    // write out data
    std::ofstream fout(outfile.c_str(), std::ios::out | std::ios::binary);
    const double T0 = 2450000.;
    std::vector<std::string> original;
    float *data   = new float [trail.npix()*trail.nspec()];
    float *errors = new float [trail.npix()*trail.nspec()];
    trail.get_data(data);
    trail.get_error(errors);
    const int nbytes = 2*sizeof(float)*trail.npix();
    for(size_t ns=0; ns<trail.nspec(); ns++){
      head.set("Record", new Subs::Hint(ns));
      head.set("HJD",    new Subs::Hdouble(T0+trail.time()[ns]));

      // Write header
      Colly::write_molly_head(fout, head, original, false);

      // Write data with Fortran numbers of bytes at start and end of each record
      if(!(fout.write((char *)&nbytes, sizeof(int))))
	throw std::string("Failed to write number of bytes at start of data record");
      if(!(fout.write((char *)(data+trail.npix()*ns),   sizeof(float[trail.npix()]))))
	throw std::string("Failed to write data");
      if(!(fout.write((char *)(errors+trail.npix()*ns), sizeof(float[trail.npix()]))))
	throw std::string("Failed to write errors");
      if(!(fout.write((char *)&nbytes, sizeof(int))))
	throw std::string("Failed to write number of bytes at end of data record");
    }
    fout.close();
    delete[] data;
    delete[] errors;

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












