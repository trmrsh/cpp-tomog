#include "trm_trail.h"

Trail::Trail(size_t npx, size_t nspc, float vp, double w0) : 
  vpix_(vp), wzero_(w0), tim_(nspc), exptim_(nspc), dat_(nspc,npx), err_(nspc,npx) {}

Trail::Trail(const std::string& file){ 
  read(file);
}

void Trail::get_data(float* arr) const {
  dat_.get(arr);
}

void Trail::set_data(float* arr) {
  dat_.set(arr);
}

void Trail::get_error(float* arr) const {
  err_.get(arr);
}

void Trail::set_error(float* arr) {
  err_.set(arr);
}

void Trail::write(const std::string& file) const{
  if(file == "-"){
    write(std::cout);
  }else{
    std::ofstream ostr(file.c_str(), std::ios::out | std::ios::binary);
    if(!ostr)
      throw Trail_Error("Trail::write -- failed to open " + file);
    write(ostr);
    ostr.close();
  }
}

void Trail::read(const std::string& file){
  if(file == "-"){
    read(std::cin);
  }else{
    std::ifstream istr(file.c_str(), std::ios::in | std::ios::binary);
    if(!istr)
      throw Trail_Error("Trail::read -- failed to open " + file);
    read(istr);
    istr.close();
  }
}

// write out trail

void Trail::write(std::ostream& ostr) const{
  
  int tflag = flag;
  ostr.write((char*)&tflag,sizeof(tflag));
  ostr.write((char*)&vpix_,sizeof(vpix_));
  ostr.write((char*)&wzero_,sizeof(wzero_));

  tim_.write(ostr);
  exptim_.write(ostr);
  dat_.write(ostr);
  err_.write(ostr);
}

void Trail::read(std::istream& istr){
  int tflag;
  istr.read((char*)&tflag,sizeof(tflag));
  if(tflag != flag) 
    throw Trail_Error("read(std::istream&) -- not a trail file");
  istr.read((char*)&vpix_,sizeof(vpix_));
  istr.read((char*)&wzero_,sizeof(wzero_));

  tim_.read(istr, false);
  exptim_.read(istr, false);
  dat_.read(istr);
  err_.read(istr);
}

bool match(const Trail& trl1, const Trail& trl2){
  return (trl1.npix() == trl2.npix() &&
	  trl1.nspec() == trl2.nspec());
}











