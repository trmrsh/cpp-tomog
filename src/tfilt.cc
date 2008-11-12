/*

!!begin
!!title  Filter a trail spectrum
!!author T.R.Marsh
!!created 13 September 2000
!!revised 22 June 2000
!!root   tfilt
!!index  tfilt
!!descr  Filter a trailed spectrum prior to back-projection
!!css   style.css
!!class  Inversion
!!class  Trailed spectra
!!head1  tfilt - filter a trailed spectrum prior to back-projection

!!emph{tfilt} applies the filter part of "filtered back-projection" to 
a set of spectra. This consists of a linearly rising part with frequency
plus a gaussian truncation for suppressing noise. The user need only specify
the FWHM of the gaussian noise suppression filter.

!!head2 Invocation

tfilt trail fwhm output!!break

!!head2 Arguments

!!table
!!arg{ trail  }{ name of a trailed spectrum.}
!!arg{ fwhm   }{ FWHM of noise suppression filter window. This is scaled in units of
cycles/pixel in which the Nyquist frequency is 0.5. Thus FWHM = 10 has little effect,
while FWHM = 0.1 has a stroing effect}
!!arg{ output }{ output file name. }
!!table

!!head2 Related commands

!!ref{tback.html}{tback}, !!ref{dtmem.html}{dtmem}

!!end

 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include "trm_subs.h"
#include "trm_constants.h"
#include "trm_input.h"
#include "trm_tomog.h"
#include "trm_trail.h"

int main(int argc, char *argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("trail",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("fwhm",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("output",  Subs::Input::LOCAL,  Subs::Input::PROMPT);

    std::string intrail;
    input.get_value("trail",  intrail, "trail", "trailed spectrum to filter");
    float fwhm;
    input.get_value("fwhm", fwhm, 0.5f, 0.00001f, 1000000.f, "FWHM of noise suppression filter (cycles/pixel)");
    std::string outfile;
    input.get_value("output", outfile, "output", "output file");

    // compute size of buffers needed 

    Trail trail(intrail);
    int npix   = trail.npix();
    int nspec  = trail.nspec();
    int nmin   = int(floor(pow(2,ceil(log(float(2*npix-1))/log(2.)))+0.1));
    int ntot   = 2*nmin;
   
    float *filter = new float[ntot];
    float *dbuff  = new float[ntot];
    float *data   = new float[npix*nspec];
    trail.get_data(data);

    const float PI2  = -1./Subs::sqr(Constants::PI);
    const float efac =  Subs::sqr(Constants::EFAC/fwhm)/2.;

    // Set convolving function which is FT of ABS(S) truncated at
    // Nyquist frequency and sampled onto same period as data. 
    // This has values 0.25 for zero frequency, 0 for
    // all even frequencies and -1/pi**2/k**2 for all odd k.
    // It is a real array. The maximum lag is dims[0]-1.
    // Note that the more obvious direct implementation is
    // avoided because of Rowland's work

    filter[0]   = 0.25;
    filter[1]   = 0.;
    float dummy;
    int i2;
    for(int ix = 1; ix < npix; ix++){
      if(2*(ix/2) == ix){
        dummy = 0.;
      }else{
        dummy = PI2/(ix*ix);
      }
      i2 = 2*ix;
      filter[i2]        = dummy; // Real part, +ve lag 
      filter[i2+1]      = 0.;    // Imaginary part, +ve lag
      filter[ntot-i2]   = dummy; // Real part, -ve lag 
      filter[ntot-i2+1] = 0.;    // Imaginary part, -ve lag 
    }

    // fill in remainder with zeroes
    for(int ix = 2*npix; ix < ntot-2*(npix-1); ix++)
      filter[ix] = 0.;

    // FFT of convolving function to get filter function
    Subs::fft(filter, ntot, 1);

    // Fold in Gaussian window function. The filter is real and symmetrical
    // so most of it is ignored
    float x;
    for(int ix = 1; ix <= nmin/2; ix++){
      i2 = 2*ix;
      x  = i2/float(nmin);
      filter[i2] *= exp(-efac*x*x);
    }

    // Now loop through spectra

    int off;
    float fac;
    
    for(int ns = 0; ns < nspec; ns++){

      // Transfer data, pad with zeroes and FFT
      
      off = ns*npix;

      for(int ix = 0; ix < npix; ix++){
	i2 = 2*ix;
	dbuff[i2]   = data[off+ix];
	dbuff[i2+1] = 0.;
      }
      for(int ix = 2*npix; ix < ntot; ix++) dbuff[ix] = 0.;

      Subs::fft(dbuff, ntot, 1);
      
      // Apply filter (which is real and symmetric) with gaussian 
      // window function, then inverse fft. First deal with zero and
      // Nyquist frequencies

      dbuff[0]      *= filter[0];
      dbuff[1]      *= filter[0];
      dbuff[nmin]   *= filter[nmin];
      dbuff[nmin+1] *= filter[nmin];

      for(int ix = 1; ix < nmin/2; ix++){
	i2    = 2*ix;
	fac   = filter[i2];
	dbuff[i2]        *= fac;
	dbuff[i2+1]      *= fac;
	dbuff[ntot-i2]   *= fac;
	dbuff[ntot-i2+1] *= fac;
      }
      Subs::fft(dbuff,ntot,-1);

      // transfer data back 
      for(int ix = 0; ix < npix; ix++)
	data[off+ix] = dbuff[2*ix]/nmin;

    }
    
    trail.set_data(data);
    trail.write(outfile);

    delete[] data;
    delete[] dbuff;
    delete[] filter;
  }

  catch(const Trail::Trail_Error& err){
    std::cerr << "Trail::Trail_Error exception: " << err << std::endl;
    exit(EXIT_FAILURE);
  }

  exit(EXIT_FAILURE);
}
