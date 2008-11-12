#ifndef  TRM_DMAP_H
#define  TRM_DMAP_H

#include <string>
#include "trm_array1d.h"
#include "trm_array2d.h"
#include "trm_buffer2d.h"

//! Represents Doppler maps
/** Dmap is able to cope with 3D Doppler maps for multiple 
 * wavelengths.
 */

class Dmap {

public:

  //! Default constructor.
  Dmap() {};

  //! Constructor of a standard Doppler map
  Dmap(int nside, float vp, float gv, double w0);

  //! Constructor of a standard Doppler map with multiple wavelengths
  Dmap(int nside, float vp, float gv, const Subs::Array1D<double>& w0);

  //! Constructor of 3D Doppler map with a single wavelengths
  Dmap(int nside, float vp, const Subs::Array1D<float>& gv, double w0);

  //! Constructor of 3D Doppler map with multiple wavelengths
  Dmap(int nside, float vp, const Subs::Array1D<float>& gv, const Subs::Array1D<double>& w0);

  //! Constructor from a named file
  Dmap(const std::string& file);
  
  //! Returns the number of pixels along a side
  int nside() const {return image_[0][0].nrow();}

  //! Returns the number of wavelengths
  int nwave() const {return image_.nrow();}

  //! Returns the number systemic velocity slices
  int ngamma() const {return image_.ncol();}

  //! Returns the total number of pixels
  int size() const {return nwave()*ngamma()*Subs::sqr(nside());}

  //! Returns the pixel size (km/s/pixel)
  float  vpix() const {return vpix_;}
  
  //! Returns the wavelengths
  Subs::Array1D<double>& wave() {return wzero_;}

  //! Returns the systemic velocities (km/s)
  Subs::Array1D<float>&  gamma() {return gamma_;}

  //! Returns a particular wavelength
  double wzero(int i) const {return wzero_[i];}

  //! Returns a particular systemic velocity (km/s)
  float gamma(int i) const {return gamma_[i];}

  //! Sets the pixel size (km/s/pixel)
  void set_vpix(float vp){vpix_=vp;}

  //! Sets a particular systemic velocity (km/s)
  void set_gamma(float gv, int i){gamma_[i]=gv;}

  //! Sets a particular systemic velocity (km/s)
  void set_wzero(double w0, int i){wzero_[i] = w0;}

  //! Returns the n-th wavelength 3D image
  Subs::Array2D<float>* operator[](int n) {return image_[n];}

  //! Returns the n-th wavelength 3D image
  const Subs::Array2D<float>* operator[](int n) const {return image_[n];}

  //! Sets the map to a constant
  Dmap& operator=(float con);

  //! Return all pixels as a single array pointer
  void get(float* arr) const;

  //! Set all pixels from a single array pointer
  void set(float* arr);

  //! Write out to a file
  void write(const std::string& file) const;

  //! Read in from a file
  void read(const std::string& file);

  //! Write out to an  opened stream
  void write(std::ostream& ostr) const;

  //! Read from an opened stream
  void read(std::istream& istr);

  //! Add a constant
  void operator+=(float con);

  //! Subtract a constant
  void operator-=(float con);

  //! Multiply by a constant
  void operator*=(float con);

  //! Divide by a constant
  void operator/=(float con);

  //! Add another image 
  void operator+=(const Dmap& dmap);

  //! Subtract another image 
  void operator-=(const Dmap& dmap);

  //! Multiply by another image 
  void operator*=(const Dmap& dmap);

  //! Divide by another image 
  void operator/=(const Dmap& dmap);

  //! Take the square root of every pixel
  void sqrt();

  //! Calculate the minimum pixel value
  float min() const;

  //! Calculate the maximum pixel value
  float max() const;

  //! Static constant to indicate file type
  const static int flag = 1235642;

  //! Error class inherited from the string class.
  class Dmap_Error : public std::string {
  public:

    //! Default constructor
    Dmap_Error() : std::string() {}

    //! Constructor from a string (e.g. an error message).
    Dmap_Error(const std::string& err) : std::string(err) {} 
  };
  

private:

  float vpix_;
  Subs::Array1D<float>  gamma_;
  Subs::Array1D<double> wzero_;
  Subs::Array2D< Subs::Array2D<float> > image_;

};

//! Calculate the maximum pixel value
float max(const Dmap& dmap);

//! Calculate the minimum pixel value
float min(const Dmap& dmap);

//! Check that formats of two maps match.
bool match(const Dmap& dmap1, const Dmap& dmap2);

//! Subtract two Dmaps
Dmap operator-(const Dmap& dmap1, const Dmap& dmap2);

//! Multiply two Dmaps
Dmap operator*(const Dmap& dmap1, const Dmap& dmap2);

//! Multiply a Dmap by a constant (constant first)
Dmap operator*(float con, const Dmap& dmap);

//! Multiply a Dmap by a constant (constant last)
Dmap operator*(const Dmap& dmap, float con);

#endif











