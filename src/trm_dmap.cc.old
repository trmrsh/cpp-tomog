#include <fstream>
#include <string>
#include "trm_subs.h"
#include "trm_dmap.h"

/** This constructs a standard 2D Doppler map of specified systemic
 * velocity and single wavelength, 
 * \param nside the number of pixels in x and y
 * \param  vp the number of km/s/pixel
 * \param  gv the gamma velocity
 * \param  w0 the rest wavelength 
 */
Dmap::Dmap(int nside, float vp, float gv, double w0) : 
  vpix_(vp), gamma_(1), wzero_(1), image_(1,1) {
  gamma_[0]  = gv;
  wzero_[0]  = w0;
  image_[0][0].resize(nside,nside);
}

/** This constructs a standard 2D Doppler map of specified systemic
 * velocity and multiple wavelengths 
 * \param nside the number of pixels in x and y
 * \param  vp the number of km/s/pixel
 * \param  gv the gamma velocity
 * \param  w0 the rest wavelength 
 */
Dmap::Dmap(int nside, float vp, float gv, const Subs::Array1D<double>& w0) : 
  vpix_(vp), gamma_(1), wzero_(w0), image_(w0.size(),1) {
  gamma_[0] = gv;
  for(int i=0; i<nwave(); i++)
    image_[i][0].resize(nside,nside);
}

// Single-wavelength, multi-gamma

Dmap::Dmap(int nside, float vp, const Subs::Array1D<float>& gv, double w0) : 
  vpix_(vp), gamma_(gv), wzero_(1), image_(1,gv.size()) {
  wzero_[0] = w0;
  for(int i=0; i<ngamma(); i++)
    image_[0][i].resize(nside,nside);
}

// Multi-wavelength, multi-gamma

Dmap::Dmap(int nside, float vp, const Subs::Array1D<float>& gv, const Subs::Array1D<double>& w0) : 
  vpix_(vp), gamma_(gv), wzero_(w0), image_(w0.size(),gv.size()) {
  for(int i=0; i<nwave(); i++)
    for(int j=0; j<ngamma(); j++)
      image_[i][j].resize(nside,nside);
}

Dmap::Dmap(const std::string& file){ 
  read(file);
}

Dmap& Dmap::operator=(float con){
  for(int i=0; i<nwave(); i++)
    for(int j=0; j<ngamma(); j++)
      image_[i][j] = con;
  return *this;
}

void Dmap::get(float* arr) const {
  for(int i=0; i<nwave(); i++)
    for(int j=0; j<ngamma(); j++){
      image_[i][j].get(arr);
      arr += image_[i][j].size();
    }
}

void Dmap::set(float* arr) {
  for(int i=0; i<nwave(); i++)
    for(int j=0; j<ngamma(); j++){
      image_[i][j].set(arr);
      arr += image_[i][j].size();
    }
}

void Dmap::write(const std::string& file) const{
  if(file == "-"){
    write(std::cout);
  }else{
    std::ofstream ostr(file.c_str(), std::ios::out | std::ios::binary);
    if(!ostr)
      throw Dmap_Error("Dmap::write -- failed to open " + file);
    write(ostr);
    ostr.close();
  }
}

void Dmap::read(const std::string& file){
  if(file == "-"){
    read(std::cin);
  }else{
    std::ifstream istr(file.c_str(), std::ios::in | std::ios::binary);
    if(!istr)
      throw Dmap_Error("Dmap::read -- failed to open " + file);
    read(istr);
    istr.close();
  }
}

// write out doppler map

void Dmap::write(std::ostream& ostr) const{

  int tflag = flag;
  ostr.write((char*)&tflag,sizeof(tflag));
  ostr.write((char*)&vpix_,sizeof(vpix_));

  wzero_.write(ostr);
  gamma_.write(ostr);

  for(int i=0; i<nwave(); i++)
    for(int j=0; j<ngamma(); j++)
      image_[i][j].write(ostr);
}

void Dmap::read(std::istream& istr){
  int tflag;
  istr.read((char*)&tflag,sizeof(tflag));
  if(tflag != flag) throw Dmap_Error("read(std::istream&) -- not a Doppler std::map file");
  istr.read((char*)&vpix_,sizeof(vpix_));

  wzero_.read(istr, false);
  gamma_.read(istr, false);

  image_.resize(wzero_.size(),gamma_.size());
  for(int i=0; i<nwave(); i++)
    for(int j=0; j<ngamma(); j++)
      image_[i][j].read(istr);

}

void Dmap::operator+=(float con){
  for(int i=0; i<nwave(); i++)
    for(int j=0; j<ngamma(); j++)
      image_[i][j] += con;
}

void Dmap::operator-=(float con){
  for(int i=0; i<nwave(); i++)
    for(int j=0; j<ngamma(); j++)
      image_[i][j] -= con;
}

void Dmap::operator*=(float con){
  for(int i=0; i<nwave(); i++)
    for(int j=0; j<ngamma(); j++)
      image_[i][j] *= con;
}

void Dmap::operator/=(float con){
  for(int i=0; i<nwave(); i++)
    for(int j=0; j<ngamma(); j++)
      image_[i][j] /= con;
}

void Dmap::operator+=(const Dmap& dmap){
  if(!match(*this,dmap)){
    throw Dmap_Error("Size mismatch in operator+=(const Dmap&)");
  }else{
    image_ += dmap.image_;
  }
}

void Dmap::operator-=(const Dmap& dmap){
  if(!match(*this,dmap)){
    throw Dmap_Error("Size mismatch in operator-=(const Dmap&)");
  }else{
    image_ -= dmap.image_;
  }
}

void Dmap::operator*=(const Dmap& dmap){
  if(!match(*this,dmap)){
    throw Dmap_Error("Size mismatch in operator*=(const Dmap&)");
  }else{
    image_ *= dmap.image_;
  }
}

void Dmap::operator/=(const Dmap& dmap){
  if(!match(*this,dmap)){
    throw Dmap_Error("Size mismatch in operator/=(const Dmap&)");
  }else{
    image_ /= dmap.image_;
  }
}

void Dmap::sqrt(){
  for(int i=0; i<nwave(); i++)
    for(int j=0; j<ngamma(); j++)
      image_[i][j].sqrt();
}

float Dmap::min() const {
  float t = image_[0][0].min(), m;
  for(int i=0; i<nwave(); i++){
    if(i){
      for(int j=0; j<ngamma(); j++)
	if(t > (m = image_[i][j].max())) t = m;
    }else{
      for(int j=1; j<ngamma(); j++)
	if(t > (m = image_[i][j].max())) t = m;
    }
  }
  return t;
}

float Dmap::max() const {
  float t = image_[0][0].max(), m;
  for(int i=0; i<nwave(); i++){
    if(i){
      for(int j=0; j<ngamma(); j++)
	if(t < (m = image_[i][j].max())) t = m;
    }else{
      for(int j=1; j<ngamma(); j++)
	if(t < (m = image_[i][j].max())) t = m;
    }
  }
  return t;
}

// non-member functions

float min(const Dmap& dmap){
  return dmap.min();
}

float max(const Dmap& dmap){
  return dmap.max();
}

bool match(const Dmap& dmap1, const Dmap& dmap2){

  if(dmap1.nwave() != dmap2.nwave() ||
     dmap1.ngamma() != dmap2.ngamma()) return false;

  for(int i=0; i<dmap1.nwave(); i++)
    for(int j=0; j<dmap1.ngamma(); j++)
      if(dmap1[i][j].nrow() != dmap2[i][j].nrow() ||
	 dmap1[i][j].ncol() != dmap2[i][j].ncol()) return false;

  return true;
}

Dmap operator-(const Dmap& dmap1, const Dmap& dmap2){
  if(!match(dmap1,dmap2)){
    throw Dmap::Dmap_Error("Size mismatch in operator-(const Dmap&, const Dmap&)");
  }else{
    Dmap dmap = dmap1;
    dmap -= dmap2;
    return dmap;
  }
}

Dmap operator*(const Dmap& dmap1, const Dmap& dmap2){
  if(!match(dmap1,dmap2)){
    throw Dmap::Dmap_Error("Size mismatch in operator*(const Dmap&, const Dmap&)");
  }else{
    Dmap dmap = dmap1;
    dmap *= dmap2;
    return dmap;
  }
}

Dmap operator*(float con, const Dmap& dmap){
  Dmap tmap = dmap;
  tmap *= con;
  return tmap;
}

Dmap operator*(const Dmap& dmap, float con){
  Dmap tmap = dmap;
  tmap *= con;
  return tmap;
}









