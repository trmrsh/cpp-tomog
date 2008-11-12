#ifndef TRM_TOMOG_H
#define TRM_TOMOG_H

#include "trm_array1d.h"

// Tomog namespace

//! Namespace of extra stuff for the tomog routines

namespace Tomog {

  //! Computes model data from a map
  void op(const float map[], const Subs::Array1D<double>& wave, 
	  const Subs::Array1D<float>& gamma, size_t nside, float vpix, 
	  float fwhm, int ndiv, int ntdiv, int npixd, int nspec, 
	  float vpixd, double waved, const Subs::Array1D<double>& time, 
	  const Subs::Array1D<float>& expose, double tzero, double period, float data[]);
  
  //! Transposed version of op
  void tr(const float data[], const Subs::Array1D<double>& wave, 
	  const Subs::Array1D<float>& gamma, size_t nside, float vpix, 
	  float fwhm, int ndiv, int ntdiv, int npixd, int nspec, 
	  float vpixd, double waved, const Subs::Array1D<double>& time, 
	  const Subs::Array1D<float>& expose, double tzero, double period, float map[]);
  
  //! Computes default image
  void gaussdef(const float input[], size_t nwave, size_t ngamma, 
		size_t nside, float fwhm, float gfwhm, float output[]);

  //! Tomog_Error is the base class for exceptions.
  class Tomog_Error : public std::string {
  public:

    //! Default constructor
    Tomog_Error() : std::string() {}

    //! Constructor from a string (e.g. an error message).
    Tomog_Error(const std::string& err) : std::string(err) {} 
  };

  //! Input_Error is an exception class for command input failures.

  /** Input_Error is a fairly heavily used class for when a command input is
   * invalid
   */
  class Input_Error : public Tomog_Error {
  public:
    //! Default constructor
    Input_Error() : Tomog_Error() {}
    //! Constructor from a string (e.g. an error message).
    Input_Error(const std::string& err) : Tomog_Error(err) {}

  };

  //! Name of environment variable which can be set to specify location of default files
  const char TOMOG_ENV[]         = "TOMOG_ENV";
  
  //! Standard name of directory for default files if environment variable not set.
  const char TOMOG_DIR[]         = ".tomog";

}

#endif
