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

//plan:
//Probably add an "astro_signal" class, more general from photometry  than can simultaneously be applied to phtometry and astrometry signals...
//Photometry should reference that...  Right now GLens is playing that role...

//Working on...
// Need to abstract the model realization from the MLdata model.
// Need a new MLmodel class which can produce lightcurve data etc
// from a vector of params.
// Then the MLdata class will produce a likelihood from a set of model data together with
// a model.  The subtlety of this is that the model may need to include some instrumental 
// parameters.
// An ideal possibility might be to allow combining models with separately defined signal
// and instrument components.  There may be a need to pull params by name,etc... need to
// think more.  Maybe there is a simpler solution...

///Base class for ml_photometry_signal
///
///This object contains information about a photometric microlensing signal model
///The signal is constructed from a GLens object, together with a Trajectory
///There are several options for controlling/modifying the form of the parameters
///applied in constructing the model, with some rough motivations.
class ML_photometry_signal : public bayes_signal{
  double r0_ref,q_ref;//2TRAJLENS:move to GLensBinary/Trajectory
  bool do_remap_r0, do_remap_q, do_log_tE;//2TRAJLENS:move to GLensBinary/Trajectory
  Trajectory *traj;
  GLens *lens;
  int idx_I0, idx_Fs, idx_q, idx_L, idx_r0, idx_phi, idx_tE, idx_tmax; //2TRAJLENS:move to GLensBinary/Traj
  stateSpace lensSpace;//2TRAJLENS:no need?
;
public:
  ML_photometry_signal(Trajectory *traj_,GLens *lens_):lens(lens_),traj(traj_){
    do_remap_r0=false;//2TRAJ:move to Trajectory
    do_remap_q=false;//2TRAJLENS:move to GLensBinary
    do_log_tE=false;//2TRAJ:move to Trajectory
    r0_ref=0;//2TRAJ:move to Trajectory
    q_ref=0;//2TRAJLENS:move to GLensBinary
    idx_I0=idx_Fs=idx_q=idx_L=idx_r0=idx_phi=idx_tE=idx_tmax=-1; //2TRAJ:clean-up
  };
  ///From bayes_signal
  //
  vector<double> get_model_signal(const state &st, const vector<double> &times)const{
    checkWorkingStateSpace();
    vector<double>params=st.get_params_vector();
    //return model_lightcurve(params,times);
    double result=0;
    double I0,Fs,q,L,r0,phi,tE,tmax;//2TRAJ
    //double I0,Fs;//2TRAJ: should look like this when we are done.
    get_model_params(params, I0,Fs,q,L,r0,phi,tE,tmax);//2TRAJ mostly eliminate this.
    vector<double> xtimes,model,modelmags;
    vector<vector<Point> > thetas;
    vector<int> indices;

    //We need to clone lens/traj before working with them so that each omp thread is working with different copies of the objects.
    GLens *worklens=clone_lens(st);
    vector<double>tEs=times;
    for(double &t : tEs)t=(t-tmax)/tE;//2TRAJ: wont need this    

    Trajectory *worktraj=clone_trajectory(worklens->getCenter(),r0,phi);
    worktraj->set_times(tEs,0);//2TRAJ: change interface to work with physical times here
    worklens->compute_trajectory(*worktraj,xtimes,thetas,indices,modelmags);

    bool burped=false;
    for(int i=0;i<indices.size();i++ ){
      double Ival = I0 - 2.5*log10(Fs*modelmags[indices[i]]+1-Fs);
      model.push_back(Ival);
      if(!isfinite(Ival)&&!burped){
	cout<<"model infinite: modelmags="<<modelmags[indices[i]]<<" I0,Fs,noise_lev,q,L,r0,phi,tE,tmax: "<<I0<<" "<<Fs<<" "<<"??????"<<" "<<q<<" "<<L<<" "<<r0<<" "<<phi<<" "<<tE<<" "<<tmax<<endl;//2TRAJ:use stateSpace for this message
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
    if(do_remap_q)idx_q=sp.requireIndex("s(1+q)");//2TRAJLENS:move to GLensBinary
    else idx_q=sp.requireIndex("logq");//2TRAJLENS:move to GLensBinary
    idx_L=sp.requireIndex("logL");//2TRAJLENS:move to GLensBinary
    if(do_remap_r0)idx_r0=sp.requireIndex("s(r0)");//2TRAJLENS:move to Trajectory
    else idx_r0=sp.requireIndex("r0");//2TRAJLENS:move to Trajectory
    idx_phi=sp.requireIndex("phi");//2TRAJLENS:move to GLensBinary
    if(do_log_tE)idx_tE=sp.requireIndex("log(tE)");//2TRAJLENS:move to Trajectory
    else idx_tE=sp.requireIndex("tE");//2TRAJLENS:move to Trajectory
    idx_tmax=sp.requireIndex("tpass");//2TRAJLENS:move to Trajectory
    haveWorkingStateSpace();
    ///Eventually want to transmit these down to constituent objects:
    ///FIXME This is a temporary version.  Should replace r0 and q rescalings with a stateSpaceTransform...
    lensSpace=lens->getObjectStateSpace();
    lens->defWorkingStateSpace(lensSpace);
  };
  
  ///Set up the output stateSpace for this object
  ///
  ///This is just an initial draft.  To be utilized in later round of development.
  stateSpace getObjectStateSpace()const{
    checkSetup();//Call this assert whenever we need options to have been processed.
    stateSpace space(9);
    string names[]={"I0","Fs","Fn","logq","logL","r0","phi","tE","tpass"};//2TRAJLENS:clean-up
    //if(use_additive_noise)names[2]="Mn";
    if(do_remap_r0)names[5]="s(r0)";//2TRAJ:move to Trajectory
    if(do_remap_q)names[3]="s(1+q)";//2TRAJLENS:move to GLensBinary
    if(do_log_tE)names[7]="log(tE)";//2TRAJ:move to Trajectory
    space.set_names(names);  
    space.set_bound(6,boundary(boundary::wrap,boundary::wrap,0,2*M_PI));//set 2-pi-wrapped space for phi.//2TRAJLENS:fix
    ///Eventually want to draw some of these up from constituent objects:
    ///eg space.add(lens.getObjectStateSpace())
    ///or space.add(transform_to_lens.inverse(lens.getObjectStateSpace()))
    return space;
  };

  void addOptions(Options &opt,const string &prefix=""){
    Optioned::addOptions(opt,prefix);
    //signal options
      addOption("log_tE","Use log10 based variable (and Gaussian prior with 1-sigma range [0:log10(tE_max)] ) for tE parameter rather than direct tE value.");//2TRAJ:move to Trajectory
      addOption("remap_r0","Use remapped r0 coordinate.");//2TRAJ:move to Trajectory
      addOption("remap_r0_ref_val","Cutoff scale for remap of r0.","2.0");//2TRAJ:move to Trajectory
      addOption("remap_q","Use remapped mass-ratio coordinate.");//2TRAJLENS:move to GLensBinary
      addOption("q0","Prior max in q (with q>1) with remapped q0. Default=1e5/","1e5");//2TRAJLENS:move to GLensBinary
  };
  void setup(){
    if(optSet("remap_r0")){//2TRAJ:move to Trajectory
      double r0_ref_val;
      *optValue("remap_r0_ref_val")>>r0_ref_val;
      remap_r0(r0_ref_val);
    }
    if(optSet("remap_q")){//2TRAJLENS:move to GLensBinary
      double q0_val;
      *optValue("q0")>>q0_val;
      remap_q(q0_val);
    }
    if(optSet("log_tE"))use_log_tE();//2TRAJ:move to Trajectory
    haveSetup();
  };    

  ///Get a new copy of the lens with appropriate parameters.
  ///This version still refers to GLensBinary parameters
  ///2TRAJLENS: This goes to GLens.
  GLens *clone_lens(const state & s)const{
    GLens *worklens;
    vector<double> params=s.get_params_vector();
    double I0,Fs,q,L,r0,phi,tE,tmax;//FIXME
    get_model_params(params, I0,Fs,q,L,r0,phi,tE,tmax);//FIXME
    worklens=lens->clone();
    state lens_state(&lensSpace,valarray<double>({q,L,0}));
    worklens->setState(lens_state);
    return worklens;
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
  //This goes to ml signal instantiation processOptions
  void remap_q(double q_ref_val=1e7){//2TRAJLENS:move to GLensBinary (or cut but preserve comment)
    do_remap_q=true;
    q_ref=q_ref_val;
  };

  ///optionally set to use logarithmic tE variable.
  //This goes to ml signal instantiation processOptions
  void use_log_tE(){
    do_log_tE=true;
  }

  ///This goes to ml instantiation of bayes_signal type
  ///2TRAJ: This is currently hard-coded with lens and trajectory parameters.  We want to be agnostic to those.
  // Here is our developing plan for addressing this:
  // Trajectory will own (tmax,tE,r0,phi, and other new params)
  // We will do conversions for rescaling these here, as needed, as with the GLens param q and Fs.
  //
  void get_model_params(const vector<double> &params, double &I0, double &Fs,double &q,double &L,double &r0,double &phi,double &tE,double &tmax)const{
    checkWorkingStateSpace();//Call this assert whenever we need the parameter index mapping.
    //Light level parameters
    //  I0 baseline unmagnitized magnitude
    //  Fs fraction of I0 light from the magnetized source
    //Earth trajectory parameterized by:
    //  time (tmax) at point of closest approach to alignment with source and lens CoM
    //  separation (r0, in Einstein ring units) of closest approach
    //  alignment angle (phi, rel to binary axis) to closest approach point
    //Intrinsic system parameters
    //  logL separation (in log10 Einstein units)
    //  q mass ratio
    I0=params[idx_I0];
    Fs=params[idx_Fs];

    //2TRAJ: wont need the rest of this func: Maybe goes to Trajectory
    //noise_lev=params[idx_noise_lev];
    double sofq=params[idx_q];//either log_q or remapped q
    double logL=params[idx_L];
    double sofr0=params[idx_r0];
    phi=params[idx_phi];
    tE=params[idx_tE];//either tE or log tE
    tmax=params[idx_tmax];
    //tmax=tmax*2;//Uncommenting this recovers pre-1-Feb-2016 def of tmax (up to mach-prec)
    
    //derived params
    if(do_log_tE)tE=pow(10.0,tE);

    //Perhaps use an alternative variable parameterizing the point of closest approach over a finite range.
    //See discussion in remap_r0() above    
    if(do_remap_r0)r0=r0_ref/sqrt(1.0/sofr0-1.0);
    else r0=sofr0;

    //2TRAJLENS: won't need the rest of this, maybe goes to GLens
    L=pow(10.0,logL);

    //Perhaps use an alternative variable parameterizing the point of closest approach over a finite range,
    //(See discussion in remap_q() above.) otherwise sofq==ln_q.
    if(do_remap_q)q=-1.0+(q_ref+1.0)/sqrt(1.0/sofq-1.0);
    else q=pow(10.0,sofq);
  };

  ///2TRAJ: Ultimately this goes into the Trajectory class as a general clone and stateSpace setting...
  ///call as traj=cone_trajectory(lens->getCenter,r0,phi)
  Trajectory *clone_trajectory(const Point &center,double r0,double phi)const{
    //compute reference position
    double vx=cos(phi), vy=sin(phi);
    double p0x=-r0*sin(phi), p0y=r0*cos(phi);
    double xcm  =  center.x;
    Trajectory *newtraj=traj->clone();
    newtraj->setup(Point(p0x+xcm,p0y), Point(vx,vy));
    return newtraj;
  };

  ///This goes to ml instantiation of bayes_signal type
  ///Change this to a standard signal function get_signal_model dep on state, data etc, returns vector.
  vector<double> old_model_lightcurve(const vector<double> &params,const vector<double>&times)const{
    double result=0;
    
    double I0,Fs,q,L,r0,phi,tE,tmax;//2TRAJ
    //double I0,Fs;//2TRAJ: should look like this when we are done.
    
    get_model_params(params, I0,Fs,q,L,r0,phi,tE,tmax);//2TRAJ
    //cout<<"model params:I0,Fs,q,L,r0,phi,tE,tmax="<<I0<<","<<Fs<<","<<q<<",\n"<<L<<","<<r0<<","<<phi<<","<<tE<<","<<tmax<<endl;
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
    //double noise_mag=I0-2.5*log10(noise_lev);
    //if(do_additive_noise)noise_mag=noise_lev; 
    vector<double> xtimes,model,modelmags;
    vector<vector<Point> > thetas;
    vector<int> indices;

    //We need to clone lens/traj before so that each omp thread is working with different copies of the objects.
    //we clone rather than copy so that we can allow derived classes.
    GLens *worklens;
    //worklens=lens->clone();
    //GLensBinary *worklens;
    //worklens=dynamic_cast<GLensBinary*>(lens->clone()); //HACK FIXME.  StateSpace based setup will obviate this cast
    //This next lines are appropriate for GLensBinary.  The generalization of this requires a state-space transformation which pulls out q and L.
    //Probably everything related to L/q should probably move into GLensBinary...
    //state lens_state(&lensSpace,valarray<double>({q,L}));//2TRAJLENS:Replace with clone_lens
    //worklens->setState(lens_state);//2TRAJLENS:Replace with clone_lens
    //worklens->setState(q,L);
    //cout<<"signal worklens:"<<worklens->print_info()<<endl;
    //worklens=new GLensBinary(q,L);
    //Trajectory traj(get_trajectory(q,L,r0,phi,tEs[0]));
    //See fixme note above.  The same issue will keep us from generalizing Trajectory, until fixed.
    //double offset=-tmax/tE;  //eliminating this changes (corrects) the meaning of t0 parameter  
    vector<double>tEs=times;
    for(double &t : tEs)t=(t-tmax)/tE;//2TRAJ: wont need this    

    Trajectory *worktraj;
    //worktraj=clone_trajectory(q,L,r0,phi,offset);
    worktraj=clone_trajectory(worklens->getCenter(),r0,phi);
    //FIXME change the following to a generic Trajectory set up based on a stateSpace setup to support more parameters.
    worktraj->set_times(tEs,0);//2TRAJ: change interface to work with physical times here
    /*
    if(debug_signal){
      int prec=cout.precision();
      cout.precision(20);
      cout<<"mlsignal Times range from "<<worktraj->t_start()<<" to "<<worktraj->t_end()<<endl;
      cout.precision(prec);
      cout<<"mlsignal Traj:"<<worktraj->print_info()<<endl;}
    */
    //cout<<"worklens="<<worklens<<endl;
    worklens->compute_trajectory(*worktraj,xtimes,thetas,indices,modelmags);

    bool burped=false;
    for(int i=0;i<indices.size();i++ ){
      double Ival = I0 - 2.5*log10(Fs*modelmags[indices[i]]+1-Fs);
      model.push_back(Ival);
      /*
      if(debug_signal){
	int prec=cout.precision();
	cout.precision(20);
	Point p=worktraj->get_obs_pos(times[i]);
	double xcm  =  (q/(1.0+q)-0.5)*L;
	cout<<i<<" "<<times[i]<<" "<<xtimes[indices[i]]<<"\nq,L="<<q<<","<<L<<" ("<<p.x<<","<<p.y<<") "<<sqrt((p.x-xcm)*(p.x-xcm)+p.y*p.y);
	cout.precision(prec);
	cout<<"\n "<<Ival<<"\n "<<modelmags[indices[i]]<<" thetas="<<endl;
	cout.precision(20);
	for(int j=0;j<thetas[indices[i]].size();j++)cout<<"  ("<<thetas[indices[i]][j].x<<","<<thetas[indices[i]][j].y<<") "<<endl;
	cout<<"]"<<endl;//debug
	cout.precision(prec);
      }
      */
      if(!isfinite(Ival)&&!burped){
	cout<<"model infinite: modelmags="<<modelmags[indices[i]]<<" I0,Fs,noise_lev,q,L,r0,phi,tE,tmax: "<<I0<<" "<<Fs<<" "<<"??????"<<" "<<q<<" "<<L<<" "<<r0<<" "<<phi<<" "<<tE<<" "<<tmax<<endl;//2TRAJ:use stateSpace for this message
	burped=true;
      }
    }
    delete worktraj, worklens;//2TRAJ: Need to check for virtual destructor.
    return model;
  };
  
public:
  //Here we always make a square window, big enough to fit the trajectory (over the specified domain) and the lens window
  void getWindow(const state &s, Point &LLcorner,Point &URcorner, double tstart=0, double tend=0, int cent=-2){
    double I0,Fs,q,L,r0,phi,tE,tmax;//2TRAJ: As in model_lightcurve
    get_model_params(s.get_params_vector(),I0,Fs,q,L,r0,phi,tE,tmax);//2TRAJ
    Point pstart(0,0),pend(0,0);
    double margin=0,width=0,x0,y0,wx,wy;
    GLens *worklens=clone_lens(s);
    if(cent>-2){//signal to use r0 as map width and center on {rminus,CoM,rplus}, when cent={0,1,2}
      x0=worklens->getCenter(cent).x;
      width=r0;
      pstart=Point(x0-width/2,-width/2);
      pend=Point(x0+width/2,+width/2);
    } else {
      double tleft=(tstart-tmax)/tE,tright=(tend-tmax)/tE;//2TRAJ:Do we use this tstart/tend?  Probably need a Trajectory function for converting physical to lens-scaled time.
      Trajectory *tr=clone_trajectory(worklens->getCenter(),r0,phi);
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
    double I0,Fs,q,L,r0,phi,tE,tmax;//2TRAJ
    get_model_params(s.get_params_vector(), I0,Fs,q,L,r0,phi,tE,tmax);//2TRAJ

    //Get square bounding box size
    cout<<"Model params: I0  Fs  q  L  r0  phi tE tmax\n"
	<<I0<<" "<<Fs<<" "<<q<<" "<<L<<" "<<r0<<" "<<phi<<" "<<tE<<" "<<tmax<<endl;//2TRAJ: Use stateSpace or cut
    //cout<<"making traj("<<q<<" "<<L<<" "<<r0<<" "<<phi<<")"<<endl;
    //cout<<" vx="<<cos(phi)<<", vy="<<sin(phi)<<endl;
    //cout<<" p0x="<<-r0*sin(phi)<<", p0y="<<r0*cos(phi)<<endl;
    //cout<<" xcm  =  "<<(q/(1.0+q)-0.5)*L<<endl;
    //double vx=cos(phi), vy=sin(phi);
    //double p0x=-r0*sin(phi), p0y=r0*cos(phi);
    ///double xcm  =  (q/(1.0+q)-0.5)*L;
    
    GLens *worklens=clone_lens(s);
    cout<<"mlsignal:dump_trajectory: worklens="<<lens->print_info()<<endl;
    double xcm  =  worklens->getCenter().x;
    Trajectory *tr=clone_trajectory(worklens->getCenter(),r0,phi);
    vector<double>tEs=times;
    for(double &t : tEs)t=(t-tmax)/tE;//2TRAJ:don't need    
    tr->set_times(tEs,0);//2TRAJ:change to physical
    cout<<"times range from "<<tr->t_start()<<" to "<<tr->t_end()<<endl;
    cout<<tr->print_info()<<endl;
    out<<"#"<<s.get_string()<<endl;
    out<<"#1.t   2. t_rel  3.x   4.y "<<endl;
    for(auto t:tEs){//2TRAJ: this will have to be physical times
      Point p=tr->get_obs_pos(t);//converted to lens-scaled time for this?
      //out<<t+tref<<" "<<t<<" "<<p.x-xcm<<" "<<p.y<<endl;
      out<<t+tref<<" "<<t<<" "<<p.x<<" "<<p.y<<endl;
    }
  };

};

#endif




