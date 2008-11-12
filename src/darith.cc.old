/*

!!begin
!!title  Adds Doppler images
!!author T.R.Marsh
!!created 09 July 2000
!!revised 22 June 2003
!!root   dadd
!!index  dadd
!!descr  Adds Doppler images
!!css   style.css
!!class  Arithematic
!!class  Doppler images
!!head1  dadd - adds Doppler images.

!!emph{dadd} adds Doppler images to each other or constants.

!!head2 Invocation

dadd input1 input2 [nwave ngamma] output!!break

!!head2 Arguments

!!table
!!arg{ input1 }{ name of Doppler image file.}
!!arg{ input2 }{ name of a Doppler image file or a constant to add to input1}
!!arg{ nwave  }{ wavelength number to process, 0 for all}
!!arg{ ngamma }{ gamma number to process, 0 for all}
!!arg{ output }{ output Doppler image}
!!table

!!head2 Related commands:

!!ref{ddiv.html}{ddiv}, !!ref{dmul.html}{dmul}, !!ref{dsub.html}{dsub}, !!ref{dset.html}{dset}

!!end

!!begin
!!title  Divides Doppler images
!!author T.R.Marsh
!!created 09 July 2000
!!revised 22 June 2003
!!root   ddiv
!!index  ddiv
!!descr  Divides Doppler images
!!css   style.css
!!class  Arithematic
!!class  Doppler images
!!head1  ddiv - divides Doppler images.

!!emph{ddiv} divides Doppler images by each other or by constants.

!!head2 Invocation

ddiv input1 input2 [nwave ngamma] output!!break

!!head2 Arguments

!!table
!!arg{ input1 }{ name of Doppler image file.}
!!arg{ input2 }{ name of a Doppler image file or a constant to divide into input1}
!!arg{ nwave  }{ wavelength number to process, 0 for all}
!!arg{ ngamma }{ gamma number to process, 0 for all}
!!arg{ output }{ output Doppler image}
!!table

!!head2 Related commands:

!!ref{dadd.html}{dadd}, !!ref{dmul.html}{dmul}, !!ref{dsub.html}{dsub}, !!ref{dset.html}{dset}

!!end

!!begin
!!title  Multiplies Doppler images
!!author T.R.Marsh
!!created 09 July 2000
!!revised 22 June 2003
!!root   dmul
!!index  dmul
!!descr  Multiplies Doppler images
!!css   style.css
!!class  Arithematic
!!class  Doppler images
!!head1  dmul - multiplies Doppler images.

!!emph{dmul} multiplies Doppler images by each other or by constants.

!!head2 Invocation

!!emph{dmul} input1 input2 [nwave ngamma] output!!break

!!head2 Arguments

!!table
!!arg{ input1 }{ name of Doppler image file.}
!!arg{ input2 }{ name of a Doppler image file or a constant to multiply into input1}
!!arg{ nwave  }{ wavelength number to process, 0 for all}
!!arg{ ngamma }{ gamma number to process, 0 for all}
!!arg{ output }{ output Doppler image}
!!table

!!head2 Related commands:

!!ref{dadd.html}{dadd}, !!ref{ddiv.html}{ddiv}, !!ref{dsub.html}{dsub}, !!ref{dset.html}{dset}

!!end

!!begin
!!title  Subtracts Doppler images
!!author T.R.Marsh
!!created 09 July 2000
!!revised 22 June 2003
!!root   dsub
!!index  dsub
!!descr  Subtracts Doppler images
!!css   style.css
!!class  Arithematic
!!class  Doppler images
!!head1  dsub - subtracts Doppler images.

!!emph{dsub} subtracts Doppler images from each other or constants from Doppler images.

!!head2 Invocation

dsub input1 input2 [nwave ngamma] output!!break

!!head2 Arguments

!!table
!!arg{ input1 }{ name of Doppler image file.}
!!arg{ input2 }{ name of a Doppler image file or a constant to subtract from input1}
!!arg{ nwave  }{ wavelength number to process, 0 for all}
!!arg{ ngamma }{ gamma number to process, 0 for all}
!!arg{ output }{ output Doppler image}
!!table

!!head2 Related commands:

!!ref{dadd.html}{dadd}, !!ref{ddiv.html}{ddiv}, !!ref{dmul.html}{dmul}, !!ref{dset.html}{dset}

!!end

!!begin
!!title  Sets Doppler images to a constant
!!author T.R.Marsh
!!created 09 July 2000
!!revised 22 June 2003
!!root   dset
!!index  dset
!!descr  Sets Doppler images to a constant
!!css   style.css
!!class  Arithematic
!!class  Doppler images
!!head1  dset - sets Doppler images to a constant

!!emph{dset} adds Doppler images to each other or constants to
Doppler images.

!!head2 Invocation

dset input1 input2 [nwave ngamma] output!!break

!!head2 Arguments

!!table
!!arg{ input1 }{ name of Doppler image file.}
!!arg{ input2 }{ name of a Doppler image file or a constant to set input1 to}
!!arg{ nwave  }{ wavelength number to process, 0 for all}
!!arg{ ngamma }{ gamma number to process, 0 for all}
!!arg{ output }{ output Doppler image}
!!table

!!head2 Related commands:

!!ref{dadd.html}{dadd}, !!ref{ddiv.html}{ddiv}, !!ref{dmul.html}{dmul}, !!ref{dsub.html}{dsub}

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

    // Set comm to name following last slash,if any.
    std::string comm  = argv[0];
    size_t slash = comm.find_last_of('/');
    if(slash != std::string::npos) comm.erase(0,slash+1);

    const int NCOM = 5;
    std::string command[NCOM] = {"dadd", "dsub", "dmul", "ddiv", "dset"};

    bool recog = false;
    for(int i=0; i<NCOM && !recog; i++) recog = (comm == command[i]);

    if(!recog) throw Tomog::Input_Error(std::string("Could not recognise command = ") + comm);

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("input1",  Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("input2",  Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("nwave",   Subs::Input::LOCAL,  Subs::Input::NOPROMPT);
    input.sign_in("ngamma",  Subs::Input::LOCAL,  Subs::Input::NOPROMPT);
    input.sign_in("output",  Subs::Input::LOCAL,  Subs::Input::PROMPT);

    std::string infile1;
    input.get_value("input1",  infile1, "input1", "first input file");
    std::string infile2;
    std::string prompt = "second input file or constant to ";
    if(comm == "dadd"){
      prompt += "add";
    }else if(comm == "dsub"){
      prompt += "subtract";
    }else if(comm == "dmul"){
      prompt += "multiply by";
    }else if(comm == "ddiv"){
      prompt += "divide by";
    }else if(comm == "dset"){
      prompt += "set first image to";
    }
    input.get_value("input2",  infile2, "input2", prompt);
    Dmap out(infile1);
    int nwave;
    input.get_value("nwave", nwave, 0, 0, out.nwave(),  "which wavelength to operate on (0 for all)");
    int ngamma;
    input.get_value("ngamma", ngamma, 0, 0, out.ngamma(), "which systemic velocity to operate on (0 for all)");
    std::string outfile;
    input.get_value("output", outfile, "output", "output file");

    std::ifstream file(infile2.c_str());
    Dmap in2;
    bool fok;
    float con;
    if(file){
      in2.read(file);
      if(!match(out,in2)){
	std::cerr << "Conflicting Doppler std::map sizes\n";
	exit(EXIT_FAILURE);
      }
      fok = true;
    }else{
      std::istringstream istr(infile2 + " ");
      istr >> con;
      if(!istr)
	throw Tomog::Input_Error("Failed either to open the second input or interpret as a constant");
      fok = false;
    }

    if(comm == "dadd"){
      if(nwave && ngamma){
	if(fok){
	  out[nwave-1][ngamma-1] += in2[nwave-1][ngamma-1];
	}else{
	  out[nwave-1][ngamma-1] += con;
	}	
      }else if(nwave){
	for(int i=0; i<out.ngamma(); i++){
	  if(fok){
	    out[nwave-1][i] += in2[nwave-1][i];
	  }else{
	    out[nwave-1][i] += con;
	  }	
	}
      }else if(ngamma){
	for(int i=0; i<out.nwave(); i++){
	  if(fok){
	    out[i][ngamma-1] += in2[i][ngamma-1];
	  }else{
	    out[i][ngamma-1] += con;
	  }	
	}
      }else{
	if(fok){
	  out += in2;
	}else{
	  out += con;
	}
      }	

    }else if(comm == "ddiv"){
	
      if(nwave && ngamma){
	if(fok){
	  out[nwave-1][ngamma-1] /= in2[nwave-1][ngamma-1];
	}else{
	  out[nwave-1][ngamma-1] /= con;
	}	
      }else if(nwave){
	for(int i=0; i<out.ngamma(); i++){
	  if(fok){
	    out[nwave-1][i] /= in2[nwave-1][i];
	  }else{
	    out[nwave-1][i] /= con;
	  }	
	}
      }else if(ngamma){
	for(int i=0; i<out.nwave(); i++){
	  if(fok){
	    out[i][ngamma-1] /= in2[i][ngamma-1];
	  }else{
	    out[i][ngamma-1] /= con;
	  }	
	}
      }else{
	if(fok){
	  out /= in2;
	}else{
	  out /= con;
	}	
      }

    }else if(comm == "dmul"){
      if(nwave && ngamma){
	if(fok){
	  out[nwave-1][ngamma-1] *= in2[nwave-1][ngamma-1];
	}else{
	  out[nwave-1][ngamma-1] *= con;
	}	
      }else if(nwave){
	for(int i=0; i<out.ngamma(); i++){
	  if(fok){
	    out[nwave-1][i] *= in2[nwave-1][i];
	  }else{
	    out[nwave-1][i] *= con;
	  }	
	}
      }else if(ngamma){
	for(int i=0; i<out.nwave(); i++){
	  if(fok){
	    out[i][ngamma-1] *= in2[i][ngamma-1];
	  }else{
	    out[i][ngamma-1] *= con;
	  }	
	}
      }else{
	if(fok){
	  out *= in2;
	}else{
	  out *= con;
	}
      }
	
    }else if(comm == "dsub"){
      if(nwave && ngamma){
	if(fok){
	  out[nwave-1][ngamma-1] -= in2[nwave-1][ngamma-1];
	}else{
	  out[nwave-1][ngamma-1] -= con;
	}	
      }else if(nwave){
	for(int i=0; i<out.ngamma(); i++){
	  if(fok){
	    out[nwave-1][i] -= in2[nwave-1][i];
	  }else{
	    out[nwave-1][i] -= con;
	  }	
	}
      }else if(ngamma){
	for(int i=0; i<out.nwave(); i++){
	  if(fok){
	    out[i][ngamma-1] -= in2[i][ngamma-1];
	  }else{
	    out[i][ngamma-1] -= con;
	  }	
	}
      }else{
	if(fok){
	  out -= in2;
	}else{
	  out -= con;
	}	
      }

    }else if(comm == "dset"){
      if(nwave && ngamma){
	if(fok){
	  out[nwave-1][ngamma-1] = in2[nwave-1][ngamma-1];
	}else{
	  out[nwave-1][ngamma-1] = con;
	}	
      }else if(nwave){
	for(int i=0; i<out.ngamma(); i++){
	  if(fok){
	    out[nwave-1][i] = in2[nwave-1][i];
	  }else{
	    out[nwave-1][i] = con;
	  }	
	}
      }else if(ngamma){
	for(int i=0; i<out.nwave(); i++){
	  if(fok){
	    out[i][ngamma-1] = in2[i][ngamma-1];
	  }else{
	    out[i][ngamma-1] = con;
	  }	
	}
      }else{
	if(fok){
	  out = in2;
	}else{
	  out = con;
	}	
      }

    }
    out.write(outfile);
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


