#ifndef  TRM_TRAIL_H
#define  TRM_TRAIL_H

#include <string>
#include <iostream>
#include <fstream>
#include "trm_array1d.h"
#include "trm_array2d.h"

//! A class to represent trailed spectra
/**
 * Trail provides a class for holding trailed spectra.
 */

class Trail {

public:

  //! Default constructor
  Trail() {};

  //! General constructor
  Trail(size_t npx, size_t nspc, float vp, double w0);

  //! Constructor from a file
  Trail(const std::string& file);

  //! Returns the number of pixels/spectrum
  size_t  npix() const {return dat_.ncol();}

  //! Returns the number of spectra
  size_t  nspec() const {return dat_.nrow();}

  //! Returns the total number of pixels
  size_t  size() const {return npix()*nspec();}

  //! Returns the km/s/pixel
  float   vpix()const {return vpix_;}

  //! Returns the rest wavelength
  double  wzero()const {return wzero_;}

  //! Sets the km/s/pixel
  void set_vpix(float vp){vpix_=vp;}

  //! Sets the rest wavelength
  void set_wzero(double w0){wzero_=w0;}

  //! Returns the data
  const Subs::Array2D<float>& data() const {return dat_;}

  //! Sets the data
  Subs::Array2D<float>& data() {return dat_;}

  //! Returns the error
  const Subs::Array2D<float>& error() const {return err_;}

  //! Sets the errors
  Subs::Array2D<float>& error() {return err_;}

  //! Returns the times
  const Subs::Array1D<double>& time() const {return tim_;}

  //! Sets the times
  Subs::Array1D<double>& time() {return tim_;}

  //! Returns the exposure times
  const Subs::Array1D<float>& expose() const {return exptim_;}

  //! Sets the exposure times
  Subs::Array1D<float>& expose() {return exptim_;}

  //! Returns the data as a standard C-style array
  void get_data(float* arr) const;

  //! Sets the data as a standard C-style array
  void set_data(float* arr);

  //! Returns the errors as a standard C-style array
  void get_error(float* arr) const;

  //! Sets the errors as a standard C-style array
  void set_error(float* arr);

  //! Writes out the spectra to a file
  void write(const std::string& file) const;

  //! Reads the spectra from a file
  void read(const std::string& file);

  //! Writes out the spectra to a stream
  void write(std::ostream& ostr) const;

  //! Reads the spectra from a stream
  void read(std::istream& istr);

  class Trail_Error : public std::string {
  public:

    //! Default constructor
    Trail_Error() : std::string() {}

    //! Constructor from a string (e.g. an error message).
    Trail_Error(const std::string& err) : std::string(err) {} 
  };

  const static int flag = 1235641;

private:

  // Pixel size, km/s
  float vpix_; 
  // Rest wavelength
  double wzero_;
  // Times for each spectrum
  Subs::Array1D<double> tim_;
  // exposure times (same units as times)
  Subs::Array1D<float>  exptim_;
  // data, errors
  Subs::Array2D<float>  dat_, err_;

};

bool match(const Trail& trl1, const Trail& trl2);

#endif











