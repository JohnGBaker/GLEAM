///Likelihood model for microlensing
///
//Written by John Baker NASA-GSFC (2014-15)
#ifndef ML_LIKE
#define ML_LIKE

using namespace std;

#include "mlsignal.hh"
#include "mldata.hh"
#include <valarray>
#include <vector>
#include <iostream>
#include<functional>
//#include <iomanip>
//#include <fstream>
#include <ctime>
#include "omp.h"
#include "options.hh"
#include "bayesian.hh"

class ML_photometry_likelihood: public bayes_likelihood{
  sampleable_probability_function * prior;
  int count;
  double total_eval_time;
  double best_post;
  state best;
  //backward compatible hack
  bool do_additive_noise;
  int idx_I0;
  int idx_Fn;
  ///The following declaration creates a stateSpace transform object
  ///There are four arguments for the stateSpaceTransformND constructor. The first is the number of dimensions, then a vector of the "in" param names,
  ///then a vector of the transformed "out" param names, then a function (declared via the c++11 lambda 'closure' function notation) defining the transformation.
  stateSpaceTransformND noise_trans{2,{"I0","Fn"},{"I0","Mn"},[](vector<double>&v){return vector<double>({v[0],v[0]-2.5*log10(v[1])});}};
  int nevery;
public:
  ML_photometry_likelihood(stateSpace *sp, ML_photometry_data *data, ML_photometry_signal *signal, sampleable_probability_function *prior=nullptr):prior(prior),bayes_likelihood(sp,data,signal){
    ///Note: here, as before, we assume that the state space is passed in.  Maybe we should be able to compute it from the signal and data though.
    do_additive_noise=false;
    best=state(sp,sp->size());
    reset();
    //noise_trans = stateSpaceTransformND(2,{"I0","Mn"},{"I0","Fn"},[](vector<double>&v){return vector<double>({v[0],v[0]-2.5*log10(v[1])});});
    idx_Fn=idx_I0=-1;
    nevery=0;
  };
  void info_every(int n){nevery=n;};
  void reset(){
    best_post=-INFINITY;
    best=best.scalar_mult(0);
    count=0;
    total_eval_time=0;
  }
  state bestState(){return best;};
  double bestPost(){return best_post;};
  /*
  double getVariance(int i, double label){
    check();
    cout<<"mllike::getVar: dataVar,signalVar="<<data->getVariance(i)<<","<<signal->getVariance(label)<<endl;
    return data->getVariance(i)+signal->getVariance(label);};
  */
  double evaluate_log(state &s){
    //#pragma omp critical
    //cout<<"Evaluating likelihood for params"<<s.get_string()<<endl;
    valarray<double>params=s.get_params();
    //clock_t tstart=clock();
    double tstart=omp_get_wtime();

    double result=log_chi_squared(s);
    double post=result;
    if(prior)post+=prior->evaluate_log(s);
    //clock_t tend=clock();
    //double eval_time = (tend-tstart)/(double)CLOCKS_PER_SEC;
    double tend=omp_get_wtime();
    double eval_time = tend-tstart;
    #pragma omp critical
    {     
      total_eval_time+=eval_time;
      count++;
      if(nevery>0&&0==count%nevery)
	cout<<"eval_time = "<<eval_time<<"  result="<<result<<" mean = "<<total_eval_time/count<<endl; 
      if(post>best_post){
        best_post=post;
        best=state(s);
      }
      //cout<<"loglike="<<result<<"<="<<maxLike<<endl;   
      if(!isfinite(result)){
        cout<<"Whoa dude, loglike is NAN! What's up with that?"<<endl;
        cout<<"params="<<s.get_string()<<endl;
	result=-INFINITY;
      }
    }
    return result;
  };

  ///from stateSpaceInterface
  void defWorkingStateSpace(const stateSpace &sp){
    checkSetup();//Call this assert whenever we need options to have been processed.
    haveWorkingStateSpace();
    checkPointers();
    signal->defWorkingStateSpace(sp);
    //backward compatible hack    
    idx_I0=sp.requireIndex("I0");
    //cout<<"do_additive_noise="<<(do_additive_noise?"true":"false")<<endl;
    if(do_additive_noise){
      idx_Fn=sp.requireIndex("Mn");
      data->defWorkingStateSpace(sp);
    } else {
      idx_Fn=sp.requireIndex("Fn");
      stateSpace st=noise_trans.transform(sp);
      data->defWorkingStateSpace(st);
    }
  };

  ///Set up the output stateSpace for this object
  ///
  ///This is just an initial draft.  To be utilized in later round of development.
  stateSpace getObjectStateSpace(){
    checkSetup();//Call this assert whenever we need options to have been processed.
    checkPointers();
    stateSpace space(1);
    space=signal->getObjectStateSpace();
    cout<<"ML_photometry_likelihood::getObjectStateSpace(): This function is not yet implemented.  Need more development of stateSpace."<<endl;
    exit(1);
    /*
    space.add(getObjectStateSpace(*data));
    ///This is just a temporary hack for backward compatibility in testing.
    if(!space.relocate("Fn",2))space.relocate("Mn",2);  
    */
    return space;
  };
  void addOptions(Options &opt,const string &prefix=""){
    Optioned::addOptions(opt,prefix);   
    addOption("additive_noise","Interpret Fn->Mn as magnitude of additive noise. Fn_max is magnitude of maximum noise level (i.e. minimum noise magnitude)");
  };
  void setup(){    
    if(optSet("additive_noise"))useAdditiveNoise();
    set_like0();
    //cout<<"setup:options="<<reportOptions()<<endl;
    //cout<<"setup:do_additive_noise="<<(do_additive_noise?"true":"false")<<endl;
    haveSetup();
  };  
  state transformDataState(const state &s)const{
    //cout<<"ML_photometrylike:transfDataSt: this="<<this<<endl;
    if(!do_additive_noise){
      state st=noise_trans.transformState(s);
      return st;
    }
    return s;
  };
  state transformSignalState(const state &s)const{return s;};
  
  void getFineGrid(int & nfine, double &tfinestart, double &tfineend)const{
    checkPointers();
    nfine=data->size()*2;
    double t0,twidth,tstart,tend;
    data->getDomainLimits(tstart,tend);
    t0=data->getFocusLabel();
    double finewidth=1.5;
    tfinestart=t0-(-tstart+tend)*finewidth/2.0;
    tfineend=t0+(-tstart+tend)*finewidth/2.0; //tfine range is twice data range centered on t0
    twidth=10;//look mostly within 10 days of peak;
    //cout<<"ML_photometry_likelihood::getFineGrid: tfs="<<tfinestart<<" < ts="<<tstart<<" < t0="<<t0<<" < te="<<tend<<" < tfe="<<tfineend<<endl;
  };
    
private:
  
  ///Reparameterize Fn as the magnitude of strictly additive noise magnitude, rather than a fractional noise level.
  void useAdditiveNoise(){
    do_additive_noise=true;
  };
  
  ///This goes to ml instantiation of bayes_signal type
  ///overloads bayes_signal fn
  virtual void write(ostream &out,state &st){oldwrite(out,st);};
  virtual void writeFine(ostream &out,state &st,int nsamples=-1, double tstart=0, double tend=0){
    if(nsamples<0)getFineGrid(nsamples,tstart,tend);
    oldwrite(out,st,nsamples,tstart,tend);};
  void oldwrite(ostream &out, state&st, int nsamples=-1, double tstart=0, double tend=0){
    checkWorkingStateSpace();//Call this assert whenever we need the parameter index mapping.
    checkPointers();
    vector<double>times;
    if(nsamples<0){
      times=data->getLabels();
    } else {
      cout<<"mllike:oldWriteFine: ns,ts,te="<<nsamples<<","<<tstart<<","<<tend<<endl;
      double delta_t=(tend-tstart)/(nsamples-1.0);
      for(int i=0;i<nsamples;i++){
	double t=tstart+i*delta_t;
	times.push_back(t);
	//vector<double> loctimes(nsamples);
	//double dt=(tend-tstart)/(nsamples-1.0)/tE;
	//double t0=(tstart-tmax)/tE;
	//cout<<"dt,t0 = "<<dt<<", "<<t0<<endl;
	//for(int i=0;i<nsamples;i++)loctimes[i]=t0+dt*i;
	//traj.set_times(loctimes,0);
	/// ts = ( t - tmax )/tE
	///    = ( tstart + i*delta_t -tmax ) / tE
	///    = ( t0*tE + i*dt*tE )/tE
	///    = t0 + i*dt = loctimes[i]

      }
    }
    //double tpk=getPeakTime();
    double tpk=data->getFocusLabel();
    double time0=data->getFocusLabel(true);
    
    //backward compatibility hack
    double noise_lev=st.get_param(idx_Fn);
    double I0=st.get_param(idx_I0);
    double noise_mag=I0-2.5*log10(noise_lev);
    if(do_additive_noise)noise_mag=noise_lev;
    //endhack
    
    /*
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
    */

    I0=st.get_param(idx_I0);
    vector<double> model=signal->get_model_signal(transformSignalState(st),times);
    vector<double> dmags=data->getDeltaValues();
    vector<double> dvar=getVariances(st);

    if(nsamples<0){
      for(int i=0;i<times.size();i++){
	double S=dvar[i];
	//if(i<10)cout<<"i="<<i<<"  S="<<S<<endl;
	if(i==0)
	  out<<"#t"<<" "<<"t_vs_pk" 
	     <<" "<<"data_mag"<<" "<<"model_mag"
	     <<" "<<"data_err"<<" "<<"model_err"
	     <<endl;
	out<<times[i]+time0<<" "<<times[i]-tpk
	   <<" "<<data->getValue(i)<<" "<<model[i]
	   <<" "<<dmags[i]<<" "<<sqrt(S)
	   <<endl;
      }
    } else {
      for(int i=0;i<times.size();i++){
	//if(i<10)cout<<"i="<<i<<"  S="<<S<<endl;
	double rtS=pow(10.0,0.4*(-noise_mag+model[i]));
	double t=times[i];
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
#endif