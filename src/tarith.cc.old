/*

!!begin
!!title  Add a trail spectrum
!!author T.R.Marsh
!!created 13 Sepember 2000
!!revised 30 June 2003
!!root   tadd
!!index  tadd
!!descr  Add a constant or another trail to a trail
!!css   style.css
!!class  Arithematic
!!class  Trailed spectra
!!head1  tadd - add a constant or another trail to a trail

!!emph{tadd} carries out addition on trailed spectra.
Errors from first trail propagate unchanged.

!!head2 Invocation

tadd input1 input2 output

!!head2 Arguments

!!table
!!arg{ input1  }{ name of a trailed spectrum.}
!!arg{ input2  }{ name of another trailed spectrum or a constant.}
!!arg{ output  }{ output file name. }
!!table

!!head2 Related commands

!!ref{tdiv.html}{ddiv}, 
!!ref{tsub.html}{dsub}, !!ref{tmul.html}{dmul},
!!ref{ddiv.html}{ddiv}, !!ref{dadd.html}{dadd},
!!ref{dsub.html}{dsub}, !!ref{dmul.html}{dmul}

!!end

!!begin
!!title  Divide a trail spectrum
!!author T.R.Marsh
!!created 13 Sepember 2000
!!revised 30 June 2003
!!root   tdiv
!!index  tdiv
!!descr  Divide a trail by a constant or another trail
!!css   style.css
!!class  Arithematic
!!class  Trailed spectra
!!head1  tdiv - divide a trail by a constant or another trail

!!emph{tdiv} carries out division on trailed spectra. The 
uncertainty array is divided by the same factors
as the data.

!!head2 Invocation

tdiv input1 input2 output

!!head2 Arguments

!!table
!!arg{ input1  }{ name of a trailed spectrum.}
!!arg{ input2  }{ name of another trailed spectrum or a constant.}
!!arg{ output  }{ output file name. }
!!table

!!head2 Related commands

!!ref{tadd.html}{dadd},
!!ref{tsub.html}{dsub}, !!ref{tmul.html}{dmul},
!!ref{ddiv.html}{ddiv}, !!ref{dadd.html}{dadd},
!!ref{dsub.html}{dsub}, !!ref{dmul.html}{dmul}

!!end

!!begin
!!title  Multiply a trail spectrum
!!author T.R.Marsh
!!created 13 Sepember 2000
!!revised 30 June 2003
!!root   tmul
!!index  tmul
!!descr  Divide a trail by a constant or another trail
!!css   style.css
!!class  Arithematic
!!class  Trailed spectra
!!head1  tmul - multiply a trail by a constant or another trail

!!emph{tmul} carries out multiplication on trailed spectra. The 
uncertainty array is multiplied by the same factors as the data.

!!head2 Invocation

tmul input1 input2 output

!!head2 Arguments

!!table
!!arg{ input1  }{ name of a trailed spectrum.}
!!arg{ input2  }{ name of another trailed spectrum or a constant.}
!!arg{ output  }{ output file name. }
!!table

!!head2 Related commands

!!ref{tdiv.html}{ddiv}, !!ref{tadd.html}{dadd},
!!ref{tsub.html}{dsub},
!!ref{ddiv.html}{ddiv}, !!ref{dadd.html}{dadd},
!!ref{dsub.html}{dsub}, !!ref{dmul.html}{dmul}

!!end

!!begin
!!title  Divide a trail spectrum
!!author T.R.Marsh
!!created 13 Sepember 2000
!!revised 30 June 2003
!!root   tsub
!!index  tsub
!!descr  Subtract a constant or another trail from a trail
!!css   style.css
!!class  Arithematic
!!class  Trailed spectra
!!head1  tsub - subtract a trail by a constant from a trail

!!emph{tsub} carries out subtraction on trailed spectra. Errors from first
trail propagate unchanged.

!!head2 Invocation

tdiv input1 input2 output

!!head2 Arguments

!!table
!!arg{ input1  }{ name of a trailed spectrum.}
!!arg{ input2  }{ name of another trailed spectrum or a constant.}
!!arg{ output  }{ output file name. }
!!table

!!head2 Related commands

!!ref{tdiv.html}{ddiv}, !!ref{tadd.html}{dadd},
!!ref{tmul.html}{dmul},
!!ref{ddiv.html}{ddiv}, !!ref{dadd.html}{dadd},
!!ref{dsub.html}{dsub}, !!ref{dmul.html}{dmul}

!!end

 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include "trm_subs.h"
#include "trm_input.h"
#include "trm_tomog.h"
#include "trm_trail.h"

int main(int argc, char* argv[]){

  try{

    // Set comm to name following last slash,if any.
    std::string comm  = argv[0];
    size_t slash = comm.find_last_of('/');
    if(slash != std::string::npos) comm.erase(0,slash+1);

    const int NCOM = 4;
    std::string command[NCOM] = {"tadd", "tsub", "tmul", "tdiv"};

    bool recog = false;
    for(int i=0; i<NCOM && !recog; i++) recog = (comm == command[i]);

    if(!recog) throw Tomog::Input_Error(std::string("Could not recognise command = ") + comm);

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("input1",  Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("input2",  Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("output",  Subs::Input::LOCAL,  Subs::Input::PROMPT);

    std::string infile1;
    input.get_value("input1",  infile1, "input1", "first input file");
    std::string infile2;
    std::string prompt = "second input file or constant to ";
    if(comm == "tadd"){
      prompt += "add";
    }else if(comm == "tsub"){
      prompt += "subtract";
    }else if(comm == "tmul"){
      prompt += "multiply by";
    }else if(comm == "tdiv"){
      prompt += "divide by";
    }
    input.get_value("input2",  infile2, "input2", prompt);
    std::string outfile;
    input.get_value("output", outfile, "output", "output file");

    Trail out(infile1);
    std::ifstream file(infile2.c_str());
    Trail in2;
    bool fok;
    float con;
    if(file){
      in2.read(infile2);
      if(!match(out,in2))
	throw Tomog::Input_Error("Conflicting trail sizes");
      fok = true;
    }else{
      std::istringstream istr(infile2 + " ");
      istr >> con;
      if(!istr)
	throw Tomog::Input_Error("Failed either to open the second input or interpret as a constant");
      fok = false;
    }

    if(comm == "tadd"){
      if(fok){
	out.data()  += in2.data();
      }else{
	out.data()  += con;
      }
	
    }else if(comm == "tdiv"){
      if(fok){
	out.data()  /= in2.data();
	out.error() /= in2.data();
      }else{
	out.data()  /= con;
	out.error() /= con;
      }	

    }else if(comm == "tmul"){
      if(fok){
	out.data()  *= in2.data();
	out.error() *= in2.data();
      }else{
	out.data()  *= con;
	out.error() *= con;
      }	

    }else if(comm == "tsub"){
      if(fok){
	out.data()  -= in2.data();
      }else{
	out.data()  -= con;
      }	

    }

    out.write(outfile);
  }

  catch(const Trail::Trail_Error& err){
    std::cerr << "Trail::Trail_Error exception: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  catch(const Tomog::Input_Error& err){
    std::cerr << "Tomog::Input_Error exception: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  catch(const std::string& err){
    std::cerr << "string exception: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}


