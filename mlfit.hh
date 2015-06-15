//Gravitational microlensing fitting/data
//Written by John G Baker NASA-GSFC (2014)

#include "glens.hh"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <valarray>
using namespace std;
//#define OLD_BUGGY

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
class MLdata {
protected:
  vector<double>times,mags,dmags;
  double time0;
public:
  vector<double>getTimes()const{return times;};
  void getTimeLimits(double &start, double &end)const{
    if(times.size()==0){
      cout<<"MLdata::getTimeLimits: Cannot get limit on empty object."<<endl;
      exit(1);
    }
    start=times.front();
    end=times.back();

  };
  double getPeakTime(bool original=false)const{
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
  vector<double>getMags()const{return mags;};
  vector<double>getDeltaMags()const{return dmags;};
  ///estimate effect of white noise on magnitude noise level.
  virtual double mag_noise_var(int i, double noise_mag){
    //if(i<10){
    //cout<<"mag_noise_var: i="<<i<<"  noise_mag="<<noise_mag<<"  dmags[i]="<<dmags[i]<<"  mags[i]="<<mags[i]<<endl;
    //cout<<"  result="<<dmags[i]*dmags[i]+pow(10.0,0.8*(-noise_mag+mags[i]))<<endl;
    //}
    //note: noise_mag = I0-2.5*log10(noise_lev) [or with do_additive_noise noise_mag=noise_lev;
    //so below: log10(sigma)=0.4(-noise_mag + mag) = mag-I0+2.5*log10(noise_lev);
    static const double logfactor=2.0*log10(2.5/log(10));
#ifdef OLD_BUGGY
    return dmags[i]*dmags[i]+pow(10.0,0.8*(-noise_mag+mags[i]));
#else 
    return dmags[i]*dmags[i]+pow(10.0,logfactor+0.8*(-noise_mag+mags[i]));
#endif
  };
  int size(){return times.size();};
  ///We include this here because data must be of a specific type, and thus the information for constructing a likelihood from model data should be part of the meaning of the data.
  virtual double loglikelihood(vector<double>modelData,double Sn=0)=0;
  virtual double evaluate_log(valarray<double> &params,bool integrate=false)=0;
  virtual void write(ostream &out,vector<double>vparams,bool integrate=false, int nsamples=-1, double tstart=0, double tend=0)=0;
  ///The following functions should probably be moved out to a "model" object
  virtual void get_model_params(valarray<double> &params, double &I0, double &Fs,double &noise_lev,double &q,double &L,double &r0,double &phi,double &tE,double &tmax)=0;
  virtual Trajectory get_trajectory(double q, double L, double r0, double phi, double tstart)=0;
};

//class for OGLEII-IV DIA data
class OGLEdata : public MLdata {
  double like0;
  double r0_ref,q_ref;
  bool do_remap_r0, do_remap_q, do_additive_noise, do_log_tE;
  //OGLE-IV:Photometry data file containing 5 columns: Hel.JD, I magnitude, magnitude error, seeing estimation (in pixels - 0.26"/pixel) and sky level.
protected:
  void set_like0(){
    like0=0;
    //A largely cosmetic adjustment to yield conventional likelihood level with noise_mag=0;
    for(int i=0;i<size();i++){
      double S=dmags[i]*dmags[i];
      like0+=log(S);
    }
    like0/=-2;
    cout<<"setting like0="<<like0<<endl;
  };
public:
  OGLEdata(string &filepath){
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
    set_like0();
    do_remap_r0=false;
    do_remap_q=false;
    do_additive_noise=false;
    do_log_tE=false;
    r0_ref=0;
    q_ref=0;
    time0=getPeakTime();
#ifdef OLD_BUGGY
    time0=0;
#endif
    for(double &t : times)t-=time0;//permanently offset times from their to put the peak at 0.
    return;
  };
  //double mag_noise_var(int i, double noise_mag){return dmags[i]*dmags[i]+pow(10.0,0.8*(-noise_mag+mags[i]));};
  ///Crop out some early data.
  void cropBefore(double tstart){
    while (times.size()>0&&times[0]<tstart){
      times.erase(times.begin());
      mags.erase(mags.begin());
      dmags.erase(dmags.begin());
    }
    set_like0();
  };
  double loglikelihood(vector<double>modelData,double noise_mag=0){
    //here we assume the model data is magnitude as function of time
    //and that dmags provides a 1-sigma error size.
    double sum=0;
    double nsum=0;
    //cout<<"size="<<size()<<endl;
    //cout<<"model.size()="<<modelData.size()<<endl;
    //cout<<"noise_mag="<<noise_mag<<endl;
    for(int i=0;i<size();i++){
      double d=modelData[i]-mags[i];
      double S=mag_noise_var(i,noise_mag);
      sum+=d*d/S;
      nsum+=log(S);  //note trying to move log outside loop can lead to overflow error.
      //if(!isfinite(sum))cout<<"infinite! i="<<i<<"  model="<<modelData[i]<<"  mag="<<mags[i]<<"  d="<<d<<" S="<<S<<endl;
    }
    //cout<<" sum="<<sum<<"  nsum="<<nsum<<endl;
    sum+=nsum;
    sum/=-2;
    //cout<<" sum="<<sum<<"  like0="<<like0<<endl;
    return sum-like0;
  };
  double evaluate_log(valarray<double> &params,bool integrate=false){
    double result=0;
    
    //Light level parameters
    //  I0 baseline unmagnitized magnitude
    //  Fs fraction of I0 light from the magnetized source
    //  fractional noise intensity at I0 baseline (It is conventional to fit this ad hoc away from the lens event, but we fit for this as part of the Bayesian analysis, one can thing of this as a simple stellar variability model, with more sophisticated models possible as well, naively anyway.)
    //Earth trajectory parameterized by:
    //  time (tmax) at point of closest approach to alignment with source and lens CoM
    //  separation (p0, in Einstein ring units) of closest approach
    //  alignment angle (phi, rel to binary axis) to closest approach point
    //Intrinsic system parameters
    //  L separation (in Einstein units)
    //  q mass ratio (reparameterized)
    double I0,noise_lev;
    I0=params[0];
    noise_lev=params[2];

    vector<double>model=model_lightcurve(params,integrate);
    //derived params
    //roughly translate intensity noise to a magnitude error level
    double noise_mag=I0-2.5*log10(noise_lev);
    if(do_additive_noise)noise_mag=noise_lev;
    //cout<<"I0="<<I0<<" noise_lev="<<noise_lev<<" noise_mag="<<noise_mag<<endl;
    //cout<<"mags at data+"<<&mags[0]<<"  mags[0]="<<mags[0]<<endl;
    result=loglikelihood(model,noise_mag);
    return result;
  };

  ///Reparameterize r0 point of closest approach to a new variable which has a finite range.
  ///Some issues are that: 1) Large r0 values are kind of similar in effect, and also less likely
  ///because we a) suppose that events with non negligible magnification have been preselected,
  ///and b) we need to somehow arrange finite prior probability in the interesting region.
  ///2) For small r0 a natural prior is to expect that the cumulative probability is proportional
  ///to cross-sectional area, ie cdf(r0) ~ r0^2.  We thus want a cdf with this form near in, but
  ///with cdf->1 as r0->infinity.
  ///Appropriately this implies a maximal pdf somewhere in the middle, which should be near the 
  ///expected sensitivity limit (large r0 means small magnification).
  ///A simple such function is s=(1-cdf)=1/(1+(r0_ref/r0)^2).  The max in the pdf will be at r0=r0_ref/sqrt(3) ~ r0_ref
  ///We realize this by choosing from the unit interval for s with r0=r0_ref/sqrt(1/s-1) (thi inverse cdf).
  ///Calling this function turns this behavior on and sets r0_ref.  To understand how to scale r0_ref,
  ///think of setting a cut_off magnification level A~A0.  The pdf should be maximal near this cut-off level
  ///where, because r0 is somewhat >1 we can think only of a single lens.  Then using the expression for
  ///single lens magnification we want A0~(1-9/r0_ref^4)^{-1} or r0_ref/sqrt(3)=(1+1/A0)^0.25.
  ///If we expect typical peak magnification at factor of 2 levels, then set r0_ref ~ 3*sqrt(3)/2 ~ 2.6.  If we
  ///expect most cases with 10% magnification, then set r0_ref ~ 1.9.  Default value is r0_ref=2.0 with pdf 
  ///peak at a0 ~ 1.25.  Note that for large r0, the pdf ~ 1/r0^3 ~ (A0-1)^0.75. Thus a factor of 2 decrease in
  ///the magnification excess only means a factor of ~ 1.6 decrease in prior or about -0.23 in log-prior, hopefully
  ///not yielding strong biases even for the marginal cases.
  void remap_r0(double r0_ref_val=2.0){
    do_remap_r0=true;
    r0_ref=r0_ref_val;
  };

  ///Reparameterize mass-ratio to a new variable which has a finite range.
  ///Allowing arbitrary mass ratio and separation, essentially all stars are some kind of multiple object system with
  ///either a comparable-mass or minor leading partner.  We make the simplifying assumption that subleading partners
  ///are irrelevant.  In that case we can ask "What kind of binary is it?" about any system.  If the answer is that
  ///we can't rule out a very small/distant partner, then it is a non-binary microlensing event in the usual sense.
  ///We shouldn't need explicit model comparison as long as our prior on the mass ratio give appropriate expectation
  ///for small versus large mass partners. A plausible expectation may be that there is roughly equal mass expectation
  ///at all secondary masses out to some cutoff beyond which the expectation decays for an integrable result.
  ///In terms of mass ratio, defining q>1, then with similar reasoning, we propose a PDF which is linear in 
  ///(q+1)^{-1} up to some cutoff (default 1e7) beyond which the PDF decreases.  For this we use the same simple 
  ///function used in the r0 scaling, s=(1-CDF)=c1/(1+(q0+1)^2/(q+1)^2). To make "CDF" a CDF the c1 needs to be
  ///chosen to get 0 at q=1, yielding c1=(1+(q0+1)^2/4), though this constant is irrelevant for us.
  ///If we allow values of q<1, (which are physically equivalent to 1/q>1 values, with a change of phi) then the
  ///normalization is different, and there is some slight enhancement near q~1, but the general model is still
  ///probably a reasonable prior
  ///Assuming we are interested in mass ratios out to the Earth-Sun ratio q~3e5, we can set q0~1e7 to peak at an 
  ///uninteresting value which we would interpret as effectively single-lens.
  void remap_q(double q_ref_val=1e7){
    do_remap_q=true;
    q_ref=q_ref_val;
  };

  ///optionally set to use logarithmic tE variable.
  void use_log_tE(){
    do_log_tE=true;
  }

  ///Reparameterize Fn as the magnitude of strictly additive noise magnitude, rather than a fractional noise level.
  void useAdditiveNoise(){
    do_additive_noise=true;
  };

  void get_model_params(valarray<double> &params, double &I0, double &Fs,double &noise_lev,double &q,double &L,double &r0,double &phi,double &tE,double &tmax){
    //Light level parameters
    //  I0 baseline unmagnitized magnitude
    //  Fs fraction of I0 light from the magnetized source
    //  fractional noise intensity at I0 baseline (It is conventional to fit this ad hoc away from the lens event, but we fit for this as part of the Bayesian analysis, one can thing of this as a simple stellar variability model, with more sophisticated models possible as well, naively anyway.)
    //Earth trajectory parameterized by:
    //  time (tmax) at point of closest approach to alignment with source and lens CoM
    //  separation (r0, in Einstein ring units) of closest approach
    //  alignment angle (phi, rel to binary axis) to closest approach point
    //Intrinsic system parameters
    //  logL separation (in log10 Einstein units)
    //  q mass ratio
    I0=params[0];
    Fs=params[1];
    noise_lev=params[2];
    double sofq=params[3];//either log_q or remapped q
    double logL=params[4];
    double sofr0=params[5];
    phi=params[6];
    tE=params[7];//either tE or log tE
    tmax=params[8];

    //derived params
    if(do_log_tE)tE=pow(10.0,tE);

    //Perhaps use an alternative variable parameterizing the point of closest approach over a finite range.
    //See discussion in remap_r0() above    
    if(do_remap_r0)r0=r0_ref/sqrt(1.0/sofr0-1.0);
    else r0=sofr0;
    L=pow(10.0,logL);

    //Perhaps use an alternative variable parameterizing the point of closest approach over a finite range,
    //(See discussion in remap_q() above.) otherwise sofq==ln_q.
    if(do_remap_q)q=-1.0+(q_ref+1.0)/sqrt(1.0/sofq-1.0);
    else q=pow(10.0,sofq);
  };

  Trajectory get_trajectory(double q, double L, double r0, double phi, double tstart){
    //In C++11, I believe this should be elided and thus no slower than if explicit as in inline...

    //compute reference position
    double vx=cos(phi), vy=sin(phi);
    //cout<<"r0="<<r0<<" q="<<q<<" L="<<L<<endl;
    //double sph=vy;
    //cout<<"r+ = "<<r0-q*L/(1+q)*sph<<" r- = "<<r0+L/(1+q)*sph<<"  alpha+ = "<<Fs/(r0*(1+q)-q*L*sph)<<" alpha- = "<<Fs/(r0*(1+q)+L*sph)<<endl;
    double p0x=-r0*sin(phi), p0y=r0*cos(phi);
    double xcm  =  (q/(1.0+q)-0.5)*L;

    //cout<<"vx="<<vx<<" vy="<<vy<<endl;
    //cout<<"tEs[0]="<<tEs[0]<<endl;
    //cout<<"p0x="<<p0x<<" p0y="<<p0y<<" xcm="<<xcm<<endl;
    //cout<<" tstart="<<tstart<<endl;
    Trajectory traj(Point(p0x+vx*tstart+xcm,p0y+vy*tstart), Point(vx,vy));


    return traj;
  };

  vector<double> model_lightcurve(valarray<double> &params,bool &integrate, int nsamples=-1, double tstart=0, double tend=0){
    double result=0;
    
    double I0,Fs,noise_lev,q,L,r0,phi,tE,tmax;
    get_model_params(params, I0,Fs,noise_lev,q,L,r0,phi,tE,tmax);

    //Comment:
    //The lightcurve model is centered at the midpoint between m1 and m2
    //This seems like it would unnecessarily couple the parameters as
    //the CoM should be the crucial reference point for an unresolved
    //binary.  In the model frame, the CoM will be located at 
    // xcm  = +/- L*( m1/(m1+m2) -1/2)  = +/- L*(m1-m2)/(m1+m2)/2
    // xcm  = +/- L/2*(1-q)/(1+q)   : (+/- depends on the convention for q together with orientation of masses)
    // we assume r0,phi relate to CM frame and transform below.                        
    // we could without loss of generality assume q>=1 or q<=1, which is the same as flipping the x axis here.
    // But the resulting hard q boundary can make trouble for MCMC., so it might be better to use a reflection 
    // boundary in logq but consistency would require reflecting phi at the same time, and we haven't 
    // implemented such a double reflection yet. Or maybe it doesn't matter...

    //roughly translate intensity noise to a magnitude error level
    double noise_mag=I0-2.5*log10(noise_lev);
    if(do_additive_noise)noise_mag=noise_lev;

    vector<double>tEs=times;
    //scale for implied tE
    for(double &t : tEs)t=(t-tmax)/tE;  
    GLensBinary lens(q,L);

    //Trajectory traj(get_trajectory(q,L,r0,phi,0*tEs[0]));
    Trajectory traj(get_trajectory(q,L,r0,phi,tEs[0]));

    if(nsamples<1){//by default use the data sample times.
      traj.set_times(tEs,0);
      //cout<<"tEs0,tPeak,tPeak(orig) = "<<tEs[0]<<", "<<getPeakTime()<<", "<<getPeakTime(true)<<endl;
    } else {       //otherwise we use the specified times
      vector<double> loctimes(nsamples);
      double dt=(tend-tstart)/(nsamples-1.0)/tE;
      double t0=(tstart-tmax)/tE;
      //cout<<"dt,t0 = "<<dt<<", "<<t0<<endl;
      for(int i=0;i<nsamples;i++)loctimes[i]=t0+dt*i;
      traj.set_times(loctimes,0);
    }
    //cout<<"mlfit Times range from "<<traj.t_start()<<" to "<<traj.t_end()<<endl;
    //cout<<"mlfit Traj:"<<traj.print_info()<<endl;
    vector<double> xtimes,model,modelmags;
    vector<vector<Point> > thetas;
    vector<int> indices;
    lens.compute_trajectory(traj,xtimes,thetas,indices,modelmags,integrate);
    //cout<<"times.size()="<<times.size()<<endl;
    //cout<<"xtimes.size()="<<xtimes.size()<<endl;
    //cout<<"indices.size()="<<indices.size()<<endl;
    //construct light curve

    /*
    ofstream outfile("likelihood.out");//debug
    outfile<<"#time xtime px py r data delta model delta_model err sum nsum"<<endl;//debug
    outfile.precision(10);//debug
    double sum=0,nsum=0;//debug
    cout<<"model params: "<<I0<<" "<<Fs<<" "<<noise_lev<<" "<<q<<" "<<L<<" "<<r0<<" "<<phi<<" "<<tE<<" "<<tmax<<endl;//debug
    cout<<traj.print_info()<<endl;
    */
    bool burped=false;
    for(int i=0;i<indices.size();i++ ){
      double Ival = I0 - 2.5*log10(Fs*modelmags[indices[i]]+1-Fs);
      model.push_back(Ival);
      /*
      { //debug block
	double err=Ival-mags[i];
	double S=mag_noise_var(i,noise_mag);
	double xcm  =  (q/(1.0+q)-0.5)*L;
	sum+=err*err/S;
	nsum+=log(S);  //note trying to move log outside loop can lead to overflow error.
	Point p=traj.get_obs_pos(times[i]);
	outfile<<times[i]<<" "<<xtimes[indices[i]]<<" "<<p.x<<" "<<p.y<<" "<<sqrt((p.x-xcm)*(p.x-xcm)+p.y*p.y)<<" "<<mags[i]<<" "<<dmags[i]<<" "<<Ival<<" "<<sqrt(S)<<" "<<err<<" "<<sum<<" "<<nsum<<endl;//debug
      } //debug
      */
      if(!isfinite(Ival)&&!burped){
	cout<<"model infinite: modelmags="<<modelmags[indices[i]]<<" I0,Fs,noise_lev,q,L,r0,phi,tE,tmax: "<<I0<<" "<<Fs<<" "<<noise_lev<<" "<<q<<" "<<L<<" "<<r0<<" "<<phi<<" "<<tE<<" "<<tmax<<endl;
	burped=true;
      }
    }
    //outfile.close();//debug
    return model;
  };

  void write(ostream &out,vector<double>vparams,bool integrate=false, int nsamples=-1, double tstart=0, double tend=0){
    
    valarray<double> params(vparams.data(), vparams.size());
    double I0,noise_lev;
    I0=params[0];
    noise_lev=params[2];
    double noise_mag=I0-2.5*log10(noise_lev);
    if(do_additive_noise)noise_mag=noise_lev;
    vector<double> model=model_lightcurve(params,integrate,nsamples,tstart,tend);
    double tpk=getPeakTime();

    if(nsamples<0){
      for(int i=0;i<size();i++){
	double S=mag_noise_var(i,noise_mag);
	//if(i<10)cout<<"i="<<i<<"  S="<<S<<endl;
	if(i==0)
	  out<<"#t"<<" "<<"t_vs_pk" 
	     <<" "<<"data_mag"<<" "<<"model_mag"
	     <<" "<<"data_eee"<<" "<<"model_eee"
	     <<endl;
	out<<times[i]+time0<<" "<<times[i]-tpk
	   <<" "<<mags[i]<<" "<<model[i]
	   <<" "<<dmags[i]<<" "<<sqrt(S)
	   <<endl;
      }
    } else {
      double delta_t=(tend-tstart)/(nsamples-1);
      for(int i=0;i<nsamples;i++){
	//if(i<10)cout<<"i="<<i<<"  S="<<S<<endl;
	double rtS=pow(10.0,0.4*(-noise_mag+model[i]));
	double t=tstart+i*delta_t;
	if(i==0)
	  out<<"#t"<<" "<<"t_vs_pk" 
	     <<" "<<"model_mag"<<" "<<"model_extra_err"
	     <<endl;
	out<<t+time0<<" "<<t-tpk
	   <<" "<<model[i]<<" "<<rtS
	   <<endl;
      }
    }
  };
};





