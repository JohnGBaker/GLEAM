//Gravitational lens equation for microlensing
//Written by John G Baker NASA-GSFC (2014-2016)

#ifndef TRAJECTORY_HH
#define TRAJECTORY_HH
#include <vector>
#include <iostream>
#include <cmath>
#include <sstream>
#include <iomanip>
#include "bayesian.hh"
using namespace std;
extern bool debug;

///Code for computing observer/source trajectories in relation to
/// microlensing light curves and images.

///First we define a struct for representing 2D points
typedef struct Point {
  double x;
  double y;
  Point():Point(0,0){};
  Point(double x, double y):x(x),y(y){};
  friend Point operator+(const Point &p1, const Point &p2);
  friend Point operator-(const Point &p1, const Point &p2);
  friend Point operator*(const Point &p1, const double &val);
}Point;
//

///Next is a class for trajectories through the observer plane
///base class implements a straight-line trajectory
class Trajectory : public bayes_component {
  //Transition notes:
  //In defining Trajectory parameters everything is done nominally in a frame with
  //the source-CM to lens-CM line-of-sight as origin in the observer plane.
  //Trajectory parameters indicate where the Sun/Earth passes through that
  //plane, thus Trajectory parameters include at least r0, [phi?], and probably
  //also tpass.  An issue is that, as of now, tpass is defined in physical
  //units while everything else about the trajectory is in Einstein units.
  //For Parallax trajectories, we also need the absolute physical time
  //this all means that we probably need tE available as trajectory parameter
  //as well, which is natural tE mainly relates the rate that the sun crosses
  //the observer plane to the Einstein units.
  //However, we should consider the possibility that we are working with data from
  //several observers, such as Earth/Spitzer...  In that case, the solar parameters
  //should be common to all trajectories. Still in our construction these parameters
  //seem relevant only/mainly to the Trajectorys and nowhere else.  An approach for
  //this, similar to our plan for time-alignment of data objects, can be to have
  //one trajectory serve as the base for the others, providing these common parameters.
  //This could be set up at in the setup (where parameters are realized anyway)
  //though a scheme for setting up at construction might also be considered.
  //On other consideration, we do not want to have phi necessarily be a parameter of Trajectory
  //since this is setting the relative angle between the lens and the observer
  //plane.  The observer plane orientation can be regarded as fixed wrt the suns motion direction.
  //The lens rotation and center offset can be implemented as an explicit coordinate transformation
  //as another layer with the GLens class.  That is, GLens internals will not all Traj::get_obs_pos
  //directly, but will call GLens::get_obs_pos which will apply any offset and (null by default)
  //rotation.  Overloaded versions of this can reorient the lens for binaries, this can be time
  //dependent for orbiting binaries.
protected:
  Point p0;
  Point v0;
  double ts;
  double tf;
  double cad;
  double toff;
  bool have_times;
  vector<double> times;
  double tE,tpass;
  double r0,phi,r0_ref,tE_max;
  bool do_remap_r0, do_log_tE;
  int idx_r0, idx_tE, idx_tpass; 
  stateSpace sp;
  bayes_frame *phys_time_ref_frame;
  double phys_time0;
  bool have_phys_time_ref; 
public:
  static bool verbose;
  Trajectory(Point pos0,Point vel0,double t_end,double cadence, double ts=0):p0(pos0),v0(vel0),cad(cadence),tf(t_end),ts(ts),toff(0),have_times(false),tE(1),tpass(0){
    typestring="Trajectory";option_name="LinearTraj";option_info="Linear relative trajectory.";
  };
  Trajectory(Point pos0,Point vel0):p0(pos0),v0(vel0),cad(1),tf(1),ts(0),have_times(false),tE(1),tpass(0),have_phys_time_ref(false),phys_time0(0){
    typestring="Trajectory";option_name="LinearTraj";option_info="Linear relative trajectory.";
  };
  Trajectory():Trajectory(Point(0,0),Point(1,0)){
    r0_ref=0;
    idx_r0=idx_tE=idx_tpass=-1;
    do_remap_r0=false;
    do_log_tE=false;
  };
  virtual ~Trajectory(){};//Need virtual destructor to allow derived class objects to be deleted from pointer to base.
  virtual Trajectory* clone(){
    Trajectory*cloned=new Trajectory(*this);
    return cloned;
    
  };
  //virtual void setup(const Point &pos0, const Point &vel0){p0=pos0;v0=vel0;};
  ///If there is an externally defined reference time then use this function to specify it before calling setup()
  ///The frame should provide an offset to yield physical time (in JD), from the nominal units
  virtual void set_JD_frame(bayes_frame &frame){
    if(have_phys_time_ref){
      cout<<"Trajectory::set_JD_frame: Cannot reset reference time frame."<<endl;
      exit(1);
    }
    phys_time_ref_frame=&frame;
    have_phys_time_ref=true;
  };
  ///Set the required eval times.  (Times are "phys" times, but "phys"=frame if tE=1,tpass=0)
  virtual void set_times(vector<double> times,double toff=0){
    vector<double>tEs=times;
    for(double &t : tEs)t=get_frame_time(t);
    this->times=tEs;
    this->toff=toff;
    ts=this->times[0];
    tf=this->times.back();
    have_times=true;
  };
  ///Return start time. (Times are "frame" times, but "phys"=frame if tE=1,tpass=0)
  virtual double t_start()const {return ts;};
  ///Return end time. (Times are "frame" times, but "phys"=frame if tE=1,tpass=0)
  virtual double t_end()const {return tf;};
  virtual double tEinstein()const {return tE;};
  virtual int Nsamples()const {if(have_times)return times.size(); else return (int)((t_end()-t_start())/cad)+1;};
  ///Return frame time of ith obs. 
  virtual double get_obs_time(int ith)const {if(have_times)return times[ith]; else return ts-toff+cad*ith;};
  virtual double get_phys_time(double frame_time)const{ return frame_time*tE+tpass;}//referenced to phys_time0;
  virtual double get_frame_time(double phys_time)const{return (phys_time-tpass)/tE;}
  ///Argument takes frame time below
  virtual Point get_obs_pos(double t)const {double x=p0.x+(t-toff)*v0.x,y=p0.y+(t-toff)*v0.y;return Point(x,y);};
  virtual Point get_obs_vel(double t)const {return v0;};
  virtual string print_info()const {ostringstream s;s<<"Trajectory({"<<p0.x<<","<<p0.y<<"},{"<<v0.x<<","<<v0.y<<"})"<<endl;return s.str();};
  //For bayes_component/stateSpaceInterface
  virtual void defWorkingStateSpace(const stateSpace &sp){
    checkSetup();//Call this assert whenever we need options to have been processed.
    if(do_remap_r0)idx_r0=sp.requireIndex("s(r0)");
    else idx_r0=sp.requireIndex("r0");
    if(do_log_tE)idx_tE=sp.requireIndex("log(tE)");

    else idx_tE=sp.requireIndex("tE");
    idx_tpass=sp.requireIndex("tpass");
    haveWorkingStateSpace();
  };
  virtual void setState(const state &s){
    checkWorkingStateSpace();
    r0=s.get_param(idx_r0);
    tE=s.get_param(idx_tE);
    tpass=s.get_param(idx_tpass);
    if(do_remap_r0)r0=r0_ref/sqrt(1.0/r0-1.0);
    if(do_log_tE)tE=pow(10.0,tE);
    //setup(Point(0,r0), Point(1,0));
    p0=Point(0,r0);
    v0=Point(1,0);
  };
  ///This overload of setup is (defacto?) part of the StateSpaceInterface
  virtual void setup(){
    //We internally us time referenced to the J2000 epoch start
    const double J2000day0=2451545.0;
    if(have_phys_time_ref){
      if(phys_time_ref_frame->registered()){
	phys_time0=phys_time_ref_frame->getRef()[0]-J2000day0;
      }
      else{
	vector<double> ref(1);
	ref[0]=phys_time0+J2000day0;
	phys_time_ref_frame->setRegister(ref);
      }
    }
    if(optSet("remap_r0")){
      *optValue("remap_r0_ref_val")>>r0_ref;
      do_remap_r0=true;
    }
    if(optSet("log_tE"))do_log_tE=true;
    *optValue("tE_max")>>tE_max;
    haveSetup();
    //Set nativeSpace;
    stateSpace space(3);
    string names[]={"r0","tE","tpass"};
    if(do_remap_r0)names[0]="s(r0)";
    if(do_log_tE)names[1]="log(tE)";
    space.set_names(names);  
    nativeSpace=space;
    sp=space;
    //set nativePrior
    double twidth=300,t0=0;
    double r0s=6.0;
    if(do_remap_r0)r0s=1.0;
    const int uni=mixed_dist_product::uniform, gauss=mixed_dist_product::gaussian, pol=mixed_dist_product::polar; 
    valarray<double>    centers((initializer_list<double>){ r0s/2.0, tE_max/2,  t0     });
    valarray<double> halfwidths((initializer_list<double>){ r0s/2.0, tE_max/2,  twidth });
    valarray<int>         types((initializer_list<int>){        uni,      uni,  gauss  });
    if(do_log_tE){
      centers[1]=log10(tE_max)/2;
      halfwidths[1]=log10(tE_max)/2;
      types[1]=gauss;
    }
    setPrior(new mixed_dist_product(&sp,types,centers,halfwidths));
  };    
  ///Explanation of options:
  ///
  ///remap_r0:
  ///Reparameterize r0 point of closest approach to a new variable which has a finite range.
  ///Some issues are that: 1) Large r0 values are kind of similar in effect, and also less likely
  ///because we a) suppose that events with non negligible magnification have been preselected,
  ///and b) we need to somehow arrange finite prior probability in the interesting region.
  ///2) For small r0 a natural prior is to expect that the cumulative probability is proportional
  ///to cross-sectional area, ie cdf(r0) ~ r0^2.  We may thus want a cdf with this form near in, but
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
  ///
  ///log_tE:
  ///Set to use logarithmic tE variable parameter;
  ///
  void addOptions(Options &opt,const string &prefix=""){
    Optioned::addOptions(opt,prefix);
    addOption("log_tE","Use log10 based variable (and Gaussian prior with 1-sigma range [0:log10(tE_max)] ) for tE parameter rather than direct tE value.");
    addOption("remap_r0","Use remapped r0 coordinate.");
    addOption("remap_r0_ref_val","Cutoff scale for remap of r0.","2.0");
    addOption("tE_max","Uniform prior max in tE. Default=100.0/","100.0");
  };

};
/////////////////////////////////////////////////////////////////////////////////////////////////////


///Here is a ParallaxTrajectory class for when the observer motion cannot be neglected.
///The nominal version of this class includes the Earth's trajectory accurate to
///1e4 km or about a few minutes motion.  The class is designed to be easily generalizable
///to other orbits by overloading the function that provides the orbital trajectory.
///
///This implementation will not be ready to use until we generalize the handling of Trajectory
///objects and their parameters overall.  Right now these are hard coded in mlsignal.hh and in
///gleam.cc.  Supporting new parameters will require eliminating most or all explicit parameter
///references in gleam.cc and thus requires alternative handling of the priors. We also need
///to do a re-think of the parameterization of the trajectory/lens parameters.
///
///Parameter rethink:
/// Previously, the trajectory parameters were interpreted in the lens-frame, with an angle parameter
/// which specifies the rotation between the lens frame and the velocity direction of the
/// Earth-sun projected to the image plane.  Since the velocity is assumed to be constant, we
/// are free to more precisely interpret this angle as orientation at the point of closest
/// approach of the Sun's path to the COM line-of-sight origin.
/// To complete the specification we need another angle relating the Sun's path to the physical
/// frame (e.g. of the ecliptic). That additional angle is a heliocentric version of the so-called
/// parallax direction (See Gould2014[arxiv:1408.0797] for somewhat relevant discussion). To give
/// this parameter a name, based on that literature, we call it \f$ \phi_\mu \f$ (ie phimu). We also
/// need a parameter relating the scale of 1AU, the length unit for our ephemeris, to the Einstein
/// units in the image plane.  By convention (see eg Gaudi2010[arxiv[1002:0332]) this parameter is
/// called \f$ \pi_E \f$. A 1AU displacement in the observer plane corresponds to a piE displacement
/// in \beta in Einstein units.  This is related to the mass and distances by:
/// piE^2 = (1AU/Drel/thetaE)^2 = (1/4)(1AU/Drel)(1AU/M)
/// where Drel = D_s*D_l/(D_s-D_l).
class ParallaxTrajectory : public Trajectory {
  shared_ptr<const sampleable_probability_function> parentPrior;//remember parent prior so we can replace nativePrior
  shared_ptr<const sampleable_probability_function> parallaxPrior;
  stateSpace PTspace;
protected:
  int idx_logpiE;
  int idx_phimu;
  ///We need to know the approximate location of the source to compute the parallax
  double source_ra, source_dec;
  double source_lon,source_lat;
  double cos_source_lon,cos_source_lat;
  double sin_source_lon,sin_source_lat;
  double piE;
  double phimu;
  double dphi,cos_dphi,sin_dphi;
public:
  virtual ~ParallaxTrajectory(){};//Need virtual destructor to allow derived class objects to be deleted from pointer to base.
  ParallaxTrajectory(){
    typestring="Trajectory";option_name="EarthTraj";option_info="Earth parallax relative trajectory.";
  };
  /*
    ParallaxTrajectory(Point pos0, Point vel0, double t_end, double caden, double source_ra=0, double source_dec=0, double t2000off=0, double ts=0):Trajectory(pos0,vel0,t_end,caden,ts){
    typestring="Trajectory";option_name="EarthTraj";option_info="Earth parallax relative trajectory.";
    equatorial2ecliptic(source_ra,source_dec,source_lat,source_lon);
    cos_source_lon=cos(source_lon);
    sin_source_lon=sin(source_lon);
    cos_source_lat=cos(source_lat);
    sin_source_lat=sin(source_lat);
    };
  */
  virtual ParallaxTrajectory* clone(){
    return new ParallaxTrajectory(*this);
  };
  //virtual void rotate(double pieE,double phimu){};
  Point get_obs_pos(double t)const override{return Trajectory::get_obs_pos(t)+get_obs_pos_offset(t);};
  Point get_obs_vel(double t)const {return Trajectory::get_obs_vel(t)+get_obs_vel_offset(t);};
  ///The following functions are to help with computing the parallax
  ///First we need to define the observer orbital motion in barycentric coordinates
  void get_barycentric_observer(double t, double &r, double &theta, double &phi){};
protected:
  ///These are internal functions needed to realize the parallax
  ///First note that we will use time standardized to in JD since J2000
  ///Define the observer trajectory position in SSB coordinates
  virtual void get_obs_pos_ssb(double t, double &x, double &y, double &z)const;
  ///Define the observer trajectory velocity in SSB coordinates
  ///For this class we define an approx Earth-Moon barycenter trajectory but that can be overloaded.
  virtual void get_obs_vel_ssb(double t, double &x, double &y, double &z)const;
  ///Using the approx location of the source, transform from ssb in source-lens line-of-sight frame.
  virtual Point ssb_los_transform(double x, double y, double z)const;
  ///Compute observer position offset from ssb in source-lens line-of-sight frame.
  virtual Point get_obs_pos_offset(double t)const{
    double x=0,y=0,z=0;
    get_obs_pos_ssb(t,x,y,z);
    Point result=ssb_los_transform(x,y,z)*piE;
    if(verbose)
#pragma omp critical 
      {
	cout<<"ParallaxTrajectory::get_obs_pos_offset: x,y,z= "<<x<<","<<y<<","<<z<<endl;
	cout<<"ParallaxTrajectory::get_obs_pos_offset: t= "<<t<<endl;
	cout<<"ParallaxTrajectory::get_obs_pos_offset: = ("<<result.x<<","<<result.y<<")  piE="<<piE<<endl;
      }
    return result;
  };
  ///Compute observer velocity offset from ssb in source-lens line-of-sight frame.
  virtual Point get_obs_vel_offset(double t)const{
    double x,y,z;
    get_obs_vel_ssb(t,x,y,z);
    return ssb_los_transform(x,y,z)*piE;
  };
  ///approximately convert equatorial to ecliptic coordinates
  void equatorial2ecliptic(double source_ra, double source_dec, double source_lat, double source_lon){
    const double degrad=180.0/M_PI;
    const double eps=23.4372/degrad; //as of 2016.0, decreasing at 0.00013/yr
    const double ceps=cos(eps),seps=sin(eps);
    double sa=sin(source_ra/degrad),ca=cos(source_ra/degrad);
    double sd=sin(source_dec/degrad),cd=cos(source_dec/degrad);
    source_lon=atan2(sa*ceps+sd/cd*seps,ca);
    source_lat=sd*ceps-cd*seps*sa;
  };
  //For bayes_component/stateSpaceInterface
  virtual void defWorkingStateSpace(const stateSpace &sp){
    checkSetup();//Call this assert whenever we need options to have been processed.
    idx_logpiE=sp.requireIndex("logpiE");
    idx_phimu=sp.requireIndex("phimu");
    Trajectory::defWorkingStateSpace(sp);//Also inherit params linear Trajectory params
  };
  virtual void setState(const state &s){
    //cout<<"in ParallaxTrajectory::setState."<<endl;
    checkWorkingStateSpace();
    piE=pow(10.0,s.get_param(idx_logpiE));
    phimu=s.get_param(idx_phimu);
    //Check this! Not sure whether this assignment is consistent with the description of phimu in the comment above.
    dphi=phimu;
    cos_dphi=cos(dphi);
    sin_dphi=sin(dphi);
    Trajectory::setState(s);
  };
  ///This overload of setup is (defacto?) part of the StateSpaceInterface
  void setup(){
    Trajectory::setup();
    if(not ( optSet("source_ra") and optSet("source_dec") )){
      cout<<"ParallaxTrajectory::setup(): Options source_ra and source_dec must be set to something."<<endl;
      exit(1);
    }
    *optValue("source_ra")>>source_ra;
    *optValue("source_dec")>>source_dec;
    equatorial2ecliptic(source_ra,source_dec,source_lat,source_lon);
    cos_source_lon=cos(source_lon);
    sin_source_lon=sin(source_lon);
    cos_source_lat=cos(source_lat);
    sin_source_lat=sin(source_lat);
    //Augment nativeSpace;
    stateSpace space(2);
    string names[]={"logpiE","phimu"};
    space.set_bound(1,boundary(boundary::wrap,boundary::wrap,0,2*M_PI));//set 2-pi-wrapped space for phimu.
    space.set_names(names);  
    PTspace=space;
    nativeSpace.attach(space);
    //set nativePrior
    double logpimin=-4.0;//set 3-sigma range
    double logpimax=1.0;
    valarray<double>    centers((initializer_list<double>){  (logpimax+logpimin)/2.0, M_PI                        });
    valarray<double> halfwidths((initializer_list<double>){  (logpimax-logpimin)/6.0, M_PI                        });
    valarray<int>         types((initializer_list<int>){ mixed_dist_product::gaussian, mixed_dist_product::uniform });
    parallaxPrior=make_shared<mixed_dist_product>(&PTspace,types,centers,halfwidths);
    parentPrior=nativePrior;
    setPrior(new independent_dist_product(&nativeSpace,parentPrior.get(),parallaxPrior.get()));//append to prior

  };    
  void addOptions(Options &opt,const string &prefix=""){
    Trajectory::addOptions(opt,prefix);
    addOption("source_ra","Source R.A. for parallax (in degrees).");
    addOption("source_dec","Source Dec. for parallax (in degrees).");
  };
  string print_info()const override{ostringstream s;s<<"ParallaxTrajectory({r0="<<r0<<", tE="<<tE<<", ra="<<source_ra<<", dec="<<source_dec<<", piE="<<piE<<", phimu="<<phimu<<"})"<<endl;return s.str();};
};

#endif
