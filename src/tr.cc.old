#include "optr.h"
#include "trm_constants.h"

void tr(const float data[], const Array1D<double>& wave, 
	const Array1D<float>& gamma, size_t nside, float vpix, 
	float fwhm, int ndiv, int npixd, int nspec, float vpixd, 
	double waved, const Array1D<double>& time, 
	double tzero, double period, float map[]){

  
  int i, j, l, off;

  int nfine = ndiv*npixd;   // number of pixels in fine pixel buffer.
  float fine[nfine];        // fine buffer
  
  // blurr array stuff
  
  int nblurr = int(3.*ndiv*fwhm/vpixd);
  int nbtot = 2*nblurr+1;
  float blurr[nbtot], sigma = fwhm/EFAC;
  float efac = sqr(vpixd/ndiv/sigma)/2., sum, add;
  
  int k;
  for(k = -nblurr, sum = 0.; k<= nblurr; k++)
    sum += (blurr[nblurr+k] = exp(-efac*k*k));

  for(k=0; k< nbtot; k++)
    blurr[k] /= sum;


  // Loop through spectra

  float scale  = ndiv*vpix/vpixd; // scale factor map/fine
  float pxscale, pyscale;         // projected scale factors
  double phase, cosp, sinp;       // phase, cosine and sine.
  float fpcon;             // fine pixel offset

  size_t xp, yp;
  int np;
  float fpoff;

  for(size_t nwave=0, moff=0; nwave<wave.size(); nwave++){
    for(size_t ngamma=0; ngamma<gamma.size(); ngamma++){
      for(xp=0; xp<nside; xp++){
	for(yp=0; yp<nside; yp++){
	  map[moff++] = 0.;
	}
      }
    }
  }
  
  for(int ns=0, doff=0; ns<nspec; ns++){

    phase   = (time[ns]-tzero)/period;
    cosp    = cos(TWOPI*phase);
    sinp    = sin(TWOPI*phase);

    pxscale = -scale*cosp;
    pyscale =  scale*sinp;

    for(k=0; k<nfine; k++) fine[k] = 0.;

    // Transpose of blurr and bin section
    
    for(i=0; i<npixd; i++){
      add = data[doff++];
      off = ndiv*i;
      for(l=off; l<off+ndiv; l++){
	for(k = 0, j=l-nblurr; k<nbtot; k++, j++)
	  if(j >= 0 && j < nfine) fine[j] += blurr[k]*add;
      }
    }

    // Transpose of projection section

    for(size_t nwave=0, moff=0; nwave<wave.size(); nwave++){
      for(size_t ngamma=0; ngamma<gamma.size(); ngamma++){

	// Compute fine pixel offset factor. This shows where
	// to add in to the fine pixel array. C = speed of light
	// Two other factor account for the centres of the arrays
      
	fpcon = ndiv*((npixd-1)/2. + gamma[ngamma]/vpixd + 
		      C*1.e-3*(1.-waved/wave[nwave]))
	  -scale*(-cosp+sinp)*(nside-1)/2. + 0.5;
      
	for(yp=0; yp<nside; yp++, fpcon+=pyscale){
	  for(xp=0, fpoff=fpcon; xp<nside; xp++, moff++, fpoff+=pxscale){
	    np  = int(floor(fpoff));
	    if(np >= 0 && np < nfine) map[moff] += fine[np];
	  }
	}
      }      
    }
  }
}
