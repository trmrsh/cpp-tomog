/*

!!begin
!!title  MEM Doppler tomography
!!author T.R.Marsh
!!date   08 July 2000
!!root   dtmem
!!index  dtmem
!!descr  carries out MEM iterations.
!!css   style.css
!!class  Doppler images
!!class  Trailed spectra
!!class  Inversion
!!head1  dtmem - carries out MEM iterations

dtmem is the main routine for carrying Chi**2-constrained maximisation
of entropy in the C++ dtom routines. It allows the user to specify a
number of iterations and a parameter representing the difference in directions
of the Chi**2 and entropy metrics.

!!head2 Invocation

dtmem map trail niter caim rmax default (blurr) tlim fwhm ndiv tzero period output!!break

!!head2 Arguments

!!table
!!arg{map}   {file containing a Doppler map}
!!arg{trail} {file containing the spectra}
!!arg{niter} {number of iterations}
!!arg{caim}  {the reduced Chi**2 to aim for}
!!arg{rmax}  {limit on the relative change of pixel values (e.g. 0.2)}
!!arg{def}   {default type: 'u' for uniform, 'g' for gaussian}
!!arg{blurr} {fwhm of blurr in pixels for gaussian default. This applies
to the X,Y directions (i.e. the images).}
!!arg{gblurr}{fwhm of blurr in pixels along the gamma axis.}
!!arg{tlim}  {limit on "test" to terminate iterations}
!!arg{fwhm}  {fwhm of local line profile (km/s)}
!!arg{ndiv}  {over-sampling factor for projections}
!!arg{ntdiv} {number of points per exposure to simulate finite exposure length}
!!arg{tzero} {zero point of ephemeris}
!!arg{period}{period of ephemeris}
!!arg{output}{output Doppler map file}
!!table

It is possible to specify the same file on output as used for
input. It is only over-written at the end and so the program can
be terminated without corrupting the file.

!!end

*/

#include <climits>
#include <cstdlib>
#include <cfloat>
#include <string>
#include "trm_subs.h"
#include "trm_input.h"
#include "trm_tomog.h"
#include "trm_dmap.h"
#include "trm_trail.h"
#include "trm_memsys.h"

// Global variables to get through to opus and tropus
namespace Dtom {
  size_t nside;
  int ndiv, ntdiv, npixd, nspec;
  float vpix, vpixd, fwhm;
  double waved, tzero, period;
  Subs::Array1D<float> gamma, expose;
  Subs::Array1D<double> wave, time;
}

void Mem::opus(const int j, const int k){

  std::cerr << "    OPUS " << j+1 << " ---> " << k+1 << std::endl;

  Tomog::op(Mem::Gbl::st+Mem::Gbl::kb[j], Dtom::wave, Dtom::gamma, 
	      Dtom::nside, Dtom::vpix, Dtom::fwhm, Dtom::ndiv, Dtom::ntdiv, 
	      Dtom::npixd, Dtom::nspec, Dtom::vpixd, Dtom::waved, 
	      Dtom::time, Dtom::expose, Dtom::tzero, Dtom::period, 
	      Mem::Gbl::st+Mem::Gbl::kb[k]);
}

void Mem::tropus(const int k, const int j){

  std::cerr << "  TROPUS " << j+1 << " <--- " << k+1 << std::endl;
  
  Tomog::tr(Mem::Gbl::st+Mem::Gbl::kb[k], Dtom::wave, Dtom::gamma, 
	      Dtom::nside, Dtom::vpix, Dtom::fwhm, Dtom::ndiv, Dtom::ntdiv, 
	      Dtom::npixd, Dtom::nspec, Dtom::vpixd, Dtom::waved, 
	      Dtom::time, Dtom::expose, Dtom::tzero, Dtom::period, 
	      Mem::Gbl::st+Mem::Gbl::kb[j]);
 
}

int main(int argc, char* argv[]){

  try{

    // Set memsys buffer
    const int MXBUFF = 15000000;
    Mem::Gbl::st = new float[MXBUFF];

    // Construct Input object
    Subs::Input input(argc, argv, Tomog::TOMOG_ENV, Tomog::TOMOG_DIR);

    // Define inputs
    input.sign_in("map",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("trail",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("niter",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("caim",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("rmax",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("default", Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("blurr",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("gblurr",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("tlim",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("fwhm",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("ndiv",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("ntdiv",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("tzero",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("period",  Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("output",  Subs::Input::LOCAL,  Subs::Input::PROMPT);

    std::string inmap;
    input.get_value("map",   inmap,   "map",   "input Doppler map");
    Dmap  map(inmap);
    std::string intrail;
    input.get_value("trail", intrail, "trail", "input trailed spectrum");
    Trail trail(intrail);
    int niter;
    input.get_value("niter", niter, 10, 1, INT_MAX, "number of iterations");
    float caim;
    input.get_value("caim", caim, 1.f, 0.00001f, FLT_MAX, "reduced chi**2 to aim for");
    float rmax;
    input.get_value("rmax", rmax, 0.2f, 0.001f, 1.f, "maximum step size");
    char def;
    input.get_value("default", def, 'u', "uUgG", "default type [u(niform), g(aussian)]");
    def = toupper(def);
    float blurr, gblurr;
    if(def == 'G'){
      input.get_value("blurr", blurr, 10.f, 0.001f, 10000.f, "FWHM blurr in X and Y (pixels)");
      input.get_value("gblurr", gblurr, 10.f, 0.001f, 10000.f, "FWHM blurr in gamma (pixels)");
    }
    float tlim;
    input.get_value("tlim",   tlim, 0.f, 0.0001f, 1.f, "limiting value of 'test' to terminate iterations");
    input.get_value("fwhm",   Dtom::fwhm, 100.f, 0.0001f, 100000.f, "FWHM of local line profile (km/s)");
    input.get_value("ndiv",   Dtom::ndiv, 1, 1, 200, "over-sampling factor for map/data computations");
    input.get_value("ntdiv",  Dtom::ntdiv, 1, 1, 200, "number of points per spectrum to simulate finite exposure times");
    input.get_value("tzero",  Dtom::tzero, 0.,  -DBL_MAX, DBL_MAX, "zero-crossing time");
    input.get_value("period", Dtom::period, 0.1, 1.e-6, DBL_MAX, "period");
    std::string outfile;
    input.get_value("output", outfile, "map", "output Doppler map");
    
    // Create and load buffers for data and model. 
    Dtom::nside  = map.nside();
    Dtom::npixd  = trail.npix();
    Dtom::nspec  = trail.nspec();

    int ndat = trail.size();
    int nmod = map.size();

    // Generate mem buffer pointers
    Mem::memcore(MXBUFF,nmod,ndat);

    Dtom::wave   = map.wave();
    Dtom::gamma  = map.gamma();
    Dtom::vpix   = map.vpix();
    Dtom::vpixd  = trail.vpix();
    Dtom::waved  = trail.wzero();
    Dtom::time   = trail.time();
    Dtom::expose = trail.expose();

    // Transfer data to mem buffer
    map.get(Mem::Gbl::st+Mem::Gbl::kb[0]);
    trail.get_data(Mem::Gbl::st+Mem::Gbl::kb[20]);
    trail.get_error(Mem::Gbl::st+Mem::Gbl::kb[21]);

    for(int i = 0; i < nmod; i++){
      if(Mem::Gbl::st[Mem::Gbl::kb[0]+i] <= 0.){
	std::cerr << "Model point " << i << " = " <<
	  Mem::Gbl::st[Mem::Gbl::kb[0]+i] << " is <= 0." << std::endl;
	exit(EXIT_FAILURE);
      }
    }

    // note that we divide by the number of data points
    // even though many may be masked. this is to ensure
    // bootstrapping works

    for(int i = 0; i < ndat; i++){
      if(Mem::Gbl::st[Mem::Gbl::kb[21]+i] > 0.)
	Mem::Gbl::st[Mem::Gbl::kb[21]+i] = 
	  2./Subs::sqr(Mem::Gbl::st[Mem::Gbl::kb[21]+i])/ndat;
    }    

    float c, test, acc=1., cnew, s, rnew, snew, sumf;
    int mode;
    if(def == 'U'){
      mode = 10;
    }else if(def == 'G'){
      mode = 30;
    }else{
      throw Tomog::Tomog_Error("Could not understand default option");
    }
    for(int it=0; it<niter; it++){
      std::cerr << "\nIteration " << it+1 << std::endl;
      if(def == 'G'){
	std::cerr << "Computing gaussian default ..." << std::endl;
	Tomog::gaussdef(Mem::Gbl::st+Mem::Gbl::kb[0],map.nwave(),map.ngamma(),
			  Dtom::nside,blurr,gblurr,Mem::Gbl::st+Mem::Gbl::kb[19]);
      }
      Mem::memprm(mode,20,caim,rmax,1.,acc,c,test,cnew,s,rnew,snew,sumf);
      if(test < tlim && c <= caim) break;
    }

    // transfer and write out map

    map.set(Mem::Gbl::st+Mem::Gbl::kb[0]);
    map.write(outfile);

    // Clear 
    delete[] Mem::Gbl::st;

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

  catch(const std::bad_alloc&){
    std::cerr << "Memory allocation error" << std::endl;
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}




