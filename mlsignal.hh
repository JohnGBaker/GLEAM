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
  bool do_remap_r0, do_log_tE;//2TRAJ:move to Trajectory
  Trajectory *traj;
  GLens *lens;
  double r0_ref;
  int idx_I0, idx_Fs;
  int idx_r0, idx_tE, idx_tmax; 
;
public:
  ML_photometry_signal(Trajectory *traj_,GLens *lens_):lens(lens_),traj(traj_){
    idx_I0=idx_Fs=-1;
    do_remap_r0=false;//2TRAJ:move to Trajectory
    r0_ref=0;//2TRAJ:move to Trajectory
    idx_r0=idx_tE=idx_tmax=-1; //2TRAJ:clean-up
    do_log_tE=false;//2TRAJ:move to Trajectory
  };

  //Produce the signal model
  vector<double> get_model_signal(const state &st, const vector<double> &times)const{
    checkWorkingStateSpace();
    vector<double>params=st.get_params_vector();
    //return model_lightcurve(params,times);
    double result=0;
    double I0,Fs,r0,tE,tmax;//2TRAJ
    //double I0,Fs;//2TRAJ: should look like this when we are done.
    get_model_params(params, I0,Fs,r0,tE,tmax);//2TRAJ mostly eliminate this.
    vector<double> xtimes,model,modelmags;
    vector<vector<Point> > thetas;
    vector<int> indices;

    //We need to clone lens/traj before working with them so that each omp thread is working with different copies of the objects.
    GLens *worklens=lens->clone();
    worklens->setState(st);

    Trajectory *worktraj=clone_trajectory(r0);
    vector<double>tEs=times;
    for(double &t : tEs)t=(t-tmax)/tE;//2TRAJ: wont need this    

    worktraj->set_times(tEs,0);//2TRAJ: change interface to work with physical times here
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
    delete worktraj, worklens;//2TRAJ: Need to check for virtual destructor.
    return model;
  };

  ///From StateSpaceInterface (via bayes_signal)
  ///
  void defWorkingStateSpace(const stateSpace &sp){
    checkSetup();//Call this assert whenever we need options to have been processed.
    ///This is how the names are currently hard-coded.  We want to have these space components be supplied by the signal/data objects
    //string names[]={"I0","Fs","Fn","logq","logL","r0","phi","tE","tpass"};
    //if(use_additive_noise)names[2]="Mn";
    //if(use_remapped_r0)names[5]="s(r0)";
    //if(use_remapped_q)names[3]="s(1+q)";    
    //if(use_log_tE)names[7]="log(tE)";
    idx_I0=sp.requireIndex("I0");
    idx_Fs=sp.requireIndex("Fs");
    if(do_remap_r0)idx_r0=sp.requireIndex("s(r0)");//2TRAJ:move to Trajectory
    else idx_r0=sp.requireIndex("r0");//2TRAJ:move to Trajectory
    if(do_log_tE)idx_tE=sp.requireIndex("log(tE)");//2TRAJ:move to Trajectory
    else idx_tE=sp.requireIndex("tE");//2TRAJ:move to Trajectory
    idx_tmax=sp.requireIndex("tpass");//2TRAJ:move to Trajectory
    haveWorkingStateSpace();
    ///Eventually want to transmit these down to constituent objects:
    ///FIXME This is a temporary version.  Should replace r0 and q rescalings with a stateSpaceTransform...
    //lensSpace=lens->getObjectStateSpace();
    lens->defWorkingStateSpace(sp);
  };
  
  ///Set up the output stateSpace for this object
  ///
  ///This is just an initial draft.  To be utilized in later round of development.
  stateSpace getObjectStateSpace()const{
    checkSetup();//Call this assert whenever we need options to have been processed.
    stateSpace space(6);
    string names[]={"I0","Fs","Fn","r0","tE","tpass"};//2TRAJ:clean-up
    //if(use_additive_noise)names[2]="Mn";
    if(do_remap_r0)names[5]="s(r0)";//2TRAJ:move to Trajectory
    if(do_log_tE)names[7]="log(tE)";//2TRAJ:move to Trajectory
    space.set_names(names);  
    space.attach(lens->getObjectStateSpace());
    ///or space.attach(transform_to_lens.inverse(lens.getObjectStateSpace()))
    return space;
  };

  void addOptions(Options &opt,const string &prefix=""){
    Optioned::addOptions(opt,prefix);
    //signal options
      addOption("log_tE","Use log10 based variable (and Gaussian prior with 1-sigma range [0:log10(tE_max)] ) for tE parameter rather than direct tE value.");//2TRAJ:move to Trajectory
      addOption("remap_r0","Use remapped r0 coordinate.");//2TRAJ:move to Trajectory
      addOption("remap_r0_ref_val","Cutoff scale for remap of r0.","2.0");//2TRAJ:move to Trajectory
  };
  void setup(){
    if(optSet("remap_r0")){//2TRAJ:move to Trajectory
      double r0_ref_val;
      *optValue("remap_r0_ref_val")>>r0_ref_val;
      remap_r0(r0_ref_val);
    }
    if(optSet("log_tE"))use_log_tE();//2TRAJ:move to Trajectory
    haveSetup();
  };    

private:

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
  //This goes to ml signal instantiation processOptions
  void remap_r0(double r0_ref_val=2.0){//2TRAJ:move to Trajectory (or cut but preserve comment)
    do_remap_r0=true;
    r0_ref=r0_ref_val;
  };

  ///optionally set to use logarithmic tE variable.
  //This goes to ml signal instantiation processOptions
  void use_log_tE(){
    do_log_tE=true;
  }

  ///2TRAJ: This is currently hard-coded with lens and trajectory parameters.  We want to be agnostic to those.
  // Here is our developing plan for addressing this:
  // Trajectory will own (tmax,tE,r0 and other new params)
  // We will do conversions for rescaling these here, as needed, as with the GLens param q and Fs.
  //
  void get_model_params(const vector<double> &params, double &I0, double &Fs,double &r0,double &tE,double &tmax)const{
    checkWorkingStateSpace();//Call this assert whenever we need the parameter index mapping.
    //Light level parameters
    //  I0 baseline unmagnitized magnitude
    //  Fs fraction of I0 light from the magnetized source
    //Earth trajectory parameterized by:
    //  time (tmax) at point of closest approach to alignment with source and lens CoM
    //  separation (r0, in Einstein ring units) of closest approach
    I0=params[idx_I0];
    Fs=params[idx_Fs];

    //2TRAJ: wont need the rest of this func: Maybe goes to Trajectory
    //noise_lev=params[idx_noise_lev];
    double sofr0=params[idx_r0];
    tE=params[idx_tE];//either tE or log tE
    tmax=params[idx_tmax];
    //tmax=tmax*2;//Uncommenting this recovers pre-1-Feb-2016 def of tmax (up to mach-prec)
    
    //derived params
    if(do_log_tE)tE=pow(10.0,tE);

    //Perhaps use an alternative variable parameterizing the point of closest approach over a finite range.
    //See discussion in remap_r0() above    
    if(do_remap_r0)r0=r0_ref/sqrt(1.0/sofr0-1.0);
    else r0=sofr0;
  };

  ///2TRAJ: Ultimately this goes into the Trajectory class as a general clone and stateSpace setting...
  Trajectory *clone_trajectory(double r0)const{
    Trajectory *newtraj=traj->clone();
    newtraj->setup(Point(0,r0), Point(1,0));
    return newtraj;
  };

public:
  GLens *clone_lens()const{return lens->clone();};
  //Here we always make a square window, big enough to fit the trajectory (over the specified domain) and the lens window
  //Points referenced in this function refer to *lens frame* //consider shifting
  void getWindow(const state &s, Point &LLcorner,Point &URcorner, double tstart=0, double tend=0, int cent=-2){
    double I0,Fs,r0,tE,tmax;//2TRAJ: As in model_lightcurve
    get_model_params(s.get_params_vector(),I0,Fs,r0,tE,tmax);//2TRAJ
    Point pstart(0,0),pend(0,0);
    double margin=0,width=0,x0,y0,wx,wy;
    GLens *worklens=lens->clone();
    worklens->setState(s);

    //We start work in the lens frame
    if(cent>-2){//signal to use r0 as map width and center on {rminus,CoM,rplus}, when cent={0,1,2}
      x0=worklens->getCenter(cent).x;
      width=r0;
      pstart=Point(x0-width/2,-width/2);
      pend=Point(x0+width/2,+width/2);
    } else {
      double tleft=(tstart-tmax)/tE,tright=(tend-tmax)/tE;//2TRAJ:Do we use this tstart/tend?  Probably need a Trajectory function for converting physical to lens-scaled time.
      Trajectory *tr=clone_trajectory(r0);
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
    double I0,Fs,r0,tE,tmax;//2TRAJ
    get_model_params(s.get_params_vector(), I0,Fs,r0,tE,tmax);//2TRAJ

    //Get square bounding box size
    //cout<<"Model params: I0  Fs  q  L  r0  phi tE tmax\n"
    //	<<I0<<" "<<Fs<<" "<<q<<" "<<L<<" "<<r0<<" "<<phi<<" "<<tE<<" "<<tmax<<endl;//2TRAJ: Use stateSpace or cut

    GLens *worklens=lens->clone();
    worklens->setState(s);
    double xcm  =  worklens->getCenter().x;
    Trajectory *tr=clone_trajectory(r0);
    vector<double>tEs=times;
    for(double &t : tEs)t=(t-tmax)/tE;//2TRAJ:don't need    
    tr->set_times(tEs,0);//2TRAJ:change to physical
    cout<<"times range from "<<tr->t_start()<<" to "<<tr->t_end()<<endl;
    cout<<tr->print_info()<<endl;
    out<<"#"<<s.get_string()<<endl;
    out<<"#1.t   2. t_rel  3.x   4.y "<<endl;
    for(auto t:tEs){//2TRAJ: this will have to be physical times
      Point p=tr->get_obs_pos(t);//converted to lens-scaled time for this?  //Note: here p comes out in traj frame.
      //out<<t+tref<<" "<<t<<" "<<p.x-xcm<<" "<<p.y<<endl;
      out<<t+tref<<" "<<t<<" "<<p.x<<" "<<p.y<<endl;
    }
  };

};

#endif




