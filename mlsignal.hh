//Gravitational microlensing signal models
//
//Written by John G Baker NASA-GSFC (2014-2015)
#ifndef MLSIGNAL_HH
#define MLSIGNAL_HH
#include "glens.hh"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <valarray>
#include "bayesian.hh"

using namespace std;
extern bool debug_signal;


///Base class for ml_photometry_signal
///
///This object contains information about a photometric microlensing signal model
///The signal is constructed from a GLens object, together with a Trajectory
///There are several options for controlling/modifying the form of the parameters
///applied in constructing the model, with some rough motivations.
class ML_photometry_signal : public bayes_signal{
  Trajectory *traj;
  GLens *lens;
  int idx_I0, idx_Fs;
  stateSpace localSpace;
  shared_ptr<const sampleable_probability_function> localPrior;
public:
  ML_photometry_signal(Trajectory *traj_,GLens *lens_):lens(lens_),traj(traj_){
    idx_I0=idx_Fs=-1;
    localPrior=nullptr;
  };
  ~ML_photometry_signal(){};
  //Produce the signal model
  vector<double> get_model_signal(const state &st, const vector<double> &times){
    checkWorkingStateSpace();
    double result=0;
    double I0,Fs;
    get_model_params(st,I0,Fs);
    vector<double> xtimes,model,modelmags;
    vector<vector<Point> > thetas;
    vector<int> indices;

    //We need to clone lens/traj before working with them so that each omp thread is working with different copies of the objects.
    GLens *worklens=lens->clone();
    worklens->setState(st);

    Trajectory *worktraj=traj->clone();
    worktraj->setState(st);
    worktraj->set_times(times);
    worklens->compute_trajectory(*worktraj,xtimes,thetas,indices,modelmags);

    bool burped=false;
    for(int i=0;i<indices.size();i++ ){
      double Ival = I0 - 2.5*log10(Fs*modelmags[indices[i]]+1-Fs);
      model.push_back(Ival);
      if(!isfinite(Ival)&&!burped){
	cout<<"model infinite: modelmags="<<modelmags[indices[i]]<<" at state="<<st.show()<<endl;
	burped=true;
      }
    }
    delete worktraj;
    delete worklens;
    return model;
  };

  ///From StateSpaceInterface (via bayes_signal)
  ///
  void defWorkingStateSpace(const stateSpace &sp){
    checkSetup();//Call this assert whenever we need options to have been processed.
    idx_I0=sp.requireIndex("I0");
    idx_Fs=sp.requireIndex("Fs");
    haveWorkingStateSpace();
    //cout<<"signal::defWSS: about to def lens"<<endl;
    lens->defWorkingStateSpace(sp);
    traj->defWorkingStateSpace(sp);
  };
  
  void addOptions(Options &opt,const string &prefix=""){
    Optioned::addOptions(opt,prefix);
  };
  void setup(){
    haveSetup();
    ///Set up the full output stateSpace for this object
    stateSpace space(2);
    string names[]={"I0","Fs"};//2TRAJ:clean-up
    space.set_names(names);  
    localSpace=space;
    nativeSpace=localSpace;
    nativeSpace.attach(*lens->getObjectStateSpace());
    nativeSpace.attach(*traj->getObjectStateSpace());
    //set the prior
    const int uni=mixed_dist_product::uniform, gauss=mixed_dist_product::gaussian, pol=mixed_dist_product::polar; 
    valarray<double>    centers((initializer_list<double>){ 18.0,    0.5});
    valarray<double> halfwidths((initializer_list<double>){  5.0,    0.5});
    valarray<int>         types((initializer_list<int>){   gauss,    uni});
    localPrior.reset(new mixed_dist_product(&localSpace,types,centers,halfwidths));
    setPrior(new independent_dist_product(&nativeSpace,localPrior.get(),lens->getObjectPrior().get(),traj->getObjectPrior().get()));
  };    

private:

  void get_model_params(const state &st, double &I0, double &Fs)const{
    checkWorkingStateSpace();//Call this assert whenever we need the parameter index mapping.
    //Light level parameters
    //  I0 baseline unmagnitized magnitude
    //  Fs fraction of I0 light from the magnetized source
    I0=st.get_param(idx_I0);
    Fs=st.get_param(idx_Fs);
  };

public:
  GLens *clone_lens()const{return lens->clone();};
  //Here we always make a square window, big enough to fit the trajectory (over the specified domain) and the lens window
  //Points referenced in this function refer to *lens frame* //consider shifting
  void getWindow(const state &s, Point &LLcorner,Point &URcorner, double tstart=0, double tend=0){//, int cent=-2){
    double I0,Fs;//,r0,tE,tmax;//2TRAJ: As in model_lightcurve
    //get_model_params(s.get_params_vector(),I0,Fs,r0,tE,tmax);//2TRAJ
    get_model_params(s,I0,Fs);
    Point pstart(0,0),pend(0,0);
    double margin=0,width=0,x0,y0,wx,wy;
    GLens *worklens=lens->clone();
    worklens->setState(s);

    //We start work in the lens frame
    {
      Trajectory *tr=traj->clone();
      tr->setState(s);
      double tleft=tr->get_frame_time(tstart),tright=tr->get_frame_time(tend);
      pstart=tr->get_obs_pos(tleft);
      pend=tr->get_obs_pos(tright);
      cout<<"making mag-map window between that fits points: ("<<pstart.x<<","<<pstart.y<<") and ("<<pend.x<<","<<pend.y<<")"<<endl;
      margin=1;
      delete tr;
    }
    delete worklens;
    width=wx=abs(pstart.x-pend.x);
    wy=abs(pstart.y-pend.y);
    if(wy>width)width=wy;
    x0=pstart.x;
    if(pend.x<x0)x0=pend.x;
    y0=pstart.y;
    if(pend.y<y0)y0=pend.y;
    width+=margin;
    y0=y0-(width-wy)/2.0;
    x0=x0-(width-wx)/2.0;
    cout<<"x0,y0,width="<<x0<<", "<<y0<<", "<<width<<endl;
    LLcorner=Point(x0,y0);
    URcorner=Point(x0+width,y0+width);
    cout<<"returning: LL=("<<LLcorner.x<<","<<LLcorner.y<<") UR=("<<URcorner.x<<","<<URcorner.y<<")"<<endl;
  };    

  ///Dump the trajectory
  ///Probably moves to trajectory eventually.
  void dump_trajectory(ostream &out, state &s, vector<double> &times, double tref){
    double I0,Fs;
    get_model_params(s, I0,Fs);

    GLens *worklens=lens->clone();
    worklens->setState(s);
    double xcm  =  worklens->getCenter().x;
    Trajectory *tr=traj->clone();
    tr->setState(s);
    tr->set_times(times);
    cout<<"times range from "<<tr->t_start()<<" to "<<tr->t_end()<<endl;
    cout<<tr->print_info()<<endl;
    out<<"#"<<s.get_string()<<endl;
    out<<"#1.t   2. t_rel  3.x   4.y "<<endl;
    for(auto tph:times){
      double t=tr->get_frame_time(tph);
      Point p=tr->get_obs_pos(t);//Note: here p comes out in traj frame.
      out<<t+tref<<" "<<t<<" "<<p.x<<" "<<p.y<<endl;
    }
  };

};

#endif




