//Gravitational microlensing data
//
//Written by John G Baker NASA-GSFC (2014-15)
#ifndef MLDATA_HH
#define MLDATA_HH
#include "glens.hh"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <valarray>
#include "bayesian.hh"

using namespace std;
//plan:
// This is functional as is, but not extensible (and not logical).
// Need to abstract the model realization from the MLdata model.
// Need a new MLmodel class which can produce lightcurve data etc
// from a vector of params.
// Then the MLdata class will produce a likelihood from a set of model data together with
// a model.  The subtlety of this is that the model may need to include some instrumental 
// parameters.
// An ideal possibility might be to allow combining models with separately defined signal
// and instrument components.  There may be a need to pull params by name,etc... need to
// think more.  Maybe there is a simpler solution...

//base class for data
class ML_photometry_data : public bayes_data{
protected:
  vector<double>&times,&mags,&dmags;
  double time0;
  double extra_noise_mag;
  int idx_Mn;

public:
  ///We relabel the generic bayes_data names as times/mags/etc...
  ML_photometry_data():bayes_data(),times(labels),mags(values),dmags(dvalues),time0(label0){
    extra_noise_mag=0;;
    idx_Mn=-1;
  };
  int size()const{return times.size();};
  virtual void getDomainLimits(double &start, double &end)const{
    check_data();
    if(times.size()==0){
      cout<<"MLdata::getTimeLimits: Cannot get limit on empty object."<<endl;
      exit(1);
    }
    start=times.front();
    end=times.back();
  };
  //double getPeakTime(bool original=false)const{
  virtual double getFocusLabel(bool original=false)const{
    check_data();
    if(original||times.size()<1)return time0;
    //we assume monotonic time and magnitude data.
    double mpk=-INFINITY;
    int ipk=0;
    for( int i=0;i<times.size();i++){
      double m=-mags[i];//minus because peak negative magnitude
      if(mpk<m){
	mpk=m;
	ipk=i;
      }
    }
    return times[ipk];
  };
  ///Crop out some early data.
  ///
  ///Permanently remove early portion of data
  virtual void cropBefore(double tstart){
    check_data();
    while (times.size()>0&&times[0]<tstart){
      times.erase(times.begin());
      mags.erase(mags.begin());
      dmags.erase(dmags.begin());
    }
  };
  //vector<double>getMags()const{return mags;};
  ///May eliminate/change to only reference the state-dependent version (now realized with mag_noise_var)
  //vector<double>getDeltaMags()const{return dmags;};
  virtual double getVariance(int i)const{
    check_data();
    if(i<0||i>size()){
      cout<<"ML_photometry_data:getVariance: Index out of range."<<endl;
      exit(1);
    }
    double var = dmags[i]*dmags[i];
    static const double logfactor=2.0*log10(2.5/log(10));
    var +=pow(10.0,logfactor+0.8*(-extra_noise_mag+mags[i]));
    return var;
  };
  /*
  ///Do we need a write function?
  virtual void write(ostream &out,vector<double>vparams,bool integrate=false, int nsamples=-1, double tstart=0, double tend=0)=0;
  //This one is for the nascent "signal" interface.
  virtual void write(ostream &out,state &st, int nsamples=-1, double tstart=0, double tend=0){
    write(out,st.get_params_vector(),true,nsamples,tstart,tend);
  };
  */
  virtual void set_model(state &st){
    bayes_data::set_model(st);
    extra_noise_mag=st.get_param(idx_Mn);
  };
  ///from stateSpaceInterface
  virtual void defWorkingStateSpace(const stateSpace &sp){
    ///This is how the names are currently hard-coded.  We want to have these space components be supplied by the signal/data objects
    //string names[]={"I0","Fs","Mn","logq","logL","r0","phi","tE","tpass"};
    idx_Mn=sp.get_index("Mn");
  };

  ///Set up the output stateSpace for this object
  ///
  ///This is just an initial draft.  To be utilized in later round of development.
  virtual stateSpace getObjectStateSpace()const{
    stateSpace space(1);
    string names[]={"Mn"};
    space.set_names(names);  
    return space;
  };
};

//class for OGLEII-IV DIA data
class ML_OGLEdata : public ML_photometry_data {
  //OGLE-IV:Photometry data file containing 5 columns: Hel.JD, I magnitude, magnitude error, seeing estimation (in pixels - 0.26"/pixel) and sky level.
public:
  ML_OGLEdata(){};
  void setup(const string &filepath){
    ifstream file(filepath.c_str());
    if(file.good()){
      string line;
      while(getline(file,line)){
	if(line[0]=='#')continue;//skip comment lines
	double t,m,dm;
	stringstream(line)>>t>>m>>dm;
	times.push_back(t);
	mags.push_back(m);
	dmags.push_back(dm);
      }
    } else {
      if(filepath.size()>0){//empty path signifies go forward without data
	cout<<"OGLEData::OGLEData: Could not open file '"<<filepath<<"'."<<endl;
	exit(1);
      }
    }
    have_data=true;
    //time0=getPeakTime();
    time0=getFocusLabel();
    for(double &t : times)t-=time0;//permanently offset times from their to put the peak at 0.
    return;
  };


};



#endif
