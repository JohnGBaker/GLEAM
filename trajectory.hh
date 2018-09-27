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

extern bool use_old_labels;

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

///Compute (two times) the area (with orientation) of triangle P0-P1-P2 
///
///This function is useful for both area and orientation calcuations of a triangle defined by 3 vertex Points.
///To allow economizing multiplications when needed for orientation or for summed partial area, it returns twice the area,
///allowing division by 2 to be done once at the end of the calculation if needed.
///
///When used for orientation it returns a positive value if p2 is to the left of the directed line from p0 to p1
///yielding zero if the points align.
inline double getTwiceTriangleArea(const Point &p0,const Point &p1,const Point &p2){
  double dx1=p1.x-p0.x;
  double dx2=p2.x-p0.x;
  double dy1=p1.y-p0.y;
  double dy2=p2.y-p0.y;
  return dx1*dy2-dx2*dy1;
};

///Determine (signed) area of a polygon
///
///The polygon is provided as a list of vertex points.  An optional origin point p0 may also be provided.
inline double getPolygonArea(const vector<Point>verts, const Point &p0=Point(0,0)){
  double doublearea=0;
  int n=verts.size();
  for(int i=0;i<n;i++){
    Point p1=verts[i];
    Point p2=verts[(i+1)%n];
    doublearea+=getTwiceTriangleArea(p0,p1,p2);
  }
  return doublearea/2.0;
};
	
///Determine (signed) area of a polygon and its centroid*area
///
///The polygon is provided as a list of vertex points.  An optional origin point p0 may also be provided.
///The value of centroid*signed_area is returned.  This avoids danger of divide-by-zero for a trivial polygon.
inline double getPolygonAreaCoM(const vector<Point>verts, Point &CoM, const Point &p0=Point(0,0)){
  double doublearea=0;
  double x=0,y=0;
  int n=verts.size();
  for(int i=0;i<n;i++){
    Point p1=verts[i];
    Point p2=verts[(i+1)%n];
    double threexc=p0.x+p1.x+p2.x;
    double threeyc=p0.y+p1.y+p2.y;
    double dm=getTwiceTriangleArea(p0,p1,p2);
    doublearea+=dm;
    x+=dm*threexc;
    y+=dm*threeyc;    
  }
  CoM=Point(x/6.0,y/6.0);
  return doublearea/2.0;
}
	
///Determine whether a point lies inside a polygon.
///
///Computes the winding number for the closed path defined by the polygon to determine
///whether the point is inside.  The polygon is defined by a vector of Points corresponding
///to the ordered set of vertices.
///
///For cases where the polygon wraps around itself, the 
inline int pointInPolygon(const Point &p, const vector<Point>&verts){
  int nwind=0;
  for(int i=0;i<verts.size();i++){ //loop over edges from v[0] to v[i+1]
    int inext=(i+1)%verts.size();
    if(verts[i].y <= p.y){ //case: segment does not begin in upper half plane 
      if(verts[inext].y > p.y){ // and segment crosses into upper half plane
	if(getTwiceTriangleArea(verts[i],verts[i+1],p) > 0){
	  //use triangle area orientation to determine if segment winds positive (CCW) from Q4 to Q1
	  nwind++;
	}
      }
    } else {  //alternative case segment does begin in upper half plane
      if(verts[inext].y <= p.y){ // and crosses into lower plane or midline;
	if(getTwiceTriangleArea(verts[i],verts[i+1],p) < 0){
	  //use triangle area orientation to determine if segment winds negative (CC) from Q1 to Q4
	  nwind--;
	}
      }
    }
  }
  return nwind;
}
	
///Next is a class for trajectories through the observer plane
///base class implements a straight-line trajectory
class Trajectory : public bayes_component {
  ///Trajectory parameters connect times to observer angular positions relative to the source-CM to lens-CM line-of-sight
  ///as origin in the observer plane.  Trajectory parameters indicate where the Sun/Earth/Satellite passes through that
  ///plane, thus Trajectory parameters include at least:
  ///      tE     Einstein crossing time  (days)
  ///      r0     SSB impact parameter   (Einstein units)
  ///      tpass  time of SSB closest approach (absolute time in days)
  ///Where SSB is Solar-system barycenter.
  ///
  ///An issue is the time referencing of tpass.  Depending on how the data is set up, times may come naturally as relative
  ///to data peak, absolute Julian days or time since start of J2000 epoch. While the time-referencing doesn't matter for the
  ///base class (corresponging to observation at the Sun), for non-trivial derived-class parallax trajectories, we also need
  ///the absolute physical time reference to associate to ephemeris information. The bayes_frame time frame provides a way to
  ///anchor these time references.
  ///
  ///The convention is that r0, tE, tpass remain with respect to the passing of the SSB and that observer plane orientation
  ///is regarded as fixed wrt the suns motion direction. (An old parameter phi allowed for actively changing the trajectory
  ///angle, and internal support for that probably should be eliminated).  Note that lens rotation and center offset are
  ///handled by the GLens object be implemented as an explicit coordinate transformation
  ///
  ///Supporting the future functionality of working with data from several observers, such as Earth/Spitzer.  In that case,
  ///the solar parameters should be common to all trajectories.  An approach for this, similar to our plan for time-alignment
  ///of data objects, can be to have one trajectory serve as the base for the others, providing these common parameters.
  ///This could be set up at in the setup (where parameters are realized anyway) though a scheme for setting up at
  ///construction might also be considered.  TBD

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
  double tpass_ref_time;
  bool allow_neg_r0;//For non-parallax, this is not allowed
public:
  static bool verbose;
  Trajectory(Point pos0,Point vel0,double t_end,double cadence, double ts=0):p0(pos0),v0(vel0),cad(cadence),tf(t_end),ts(ts),toff(0),have_times(false),tE(1),tpass(0),tpass_ref_time(0),allow_neg_r0(false){
    typestring="Trajectory";option_name="LinearTraj";option_info="Linear relative trajectory.";
  };
  Trajectory(Point pos0,Point vel0):p0(pos0),v0(vel0),cad(1),tf(1),ts(0),have_times(false),tE(1),tpass(0),have_phys_time_ref(false),phys_time0(0),tpass_ref_time(0),allow_neg_r0(false){
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
    if(use_old_labels){
      if(do_remap_r0)idx_r0=sp.requireIndex("s(r0)");
      else idx_r0=sp.requireIndex("r0");
    } else {
      if(do_remap_r0)idx_r0=sp.requireIndex("f(u0)");
      else idx_r0=sp.requireIndex("u0");
    }
    if(do_log_tE)idx_tE=sp.requireIndex("log(tE)");

    else idx_tE=sp.requireIndex("tE");
    idx_tpass=sp.requireIndex("tpass");
    haveWorkingStateSpace();
  };
  virtual void setState(const state &s){
    checkWorkingStateSpace();
    r0=s.get_param(idx_r0);
    tE=s.get_param(idx_tE);
    tpass=s.get_param(idx_tpass)-tpass_ref_time;
    if(do_remap_r0)r0=r0_ref/sqrt(1.0/r0-1.0);
    if(do_log_tE)tE=pow(10.0,tE);
    //setup(Point(0,r0), Point(1,0));
    p0=Point(0,r0);
    v0=Point(1,0);
  };
  ///This overload of setup is (defacto?) part of the StateSpaceInterface
  virtual void setup(){
    //We internally use time referenced to the J2000 epoch start
    const double J2000day0=2451545.0;
    if(have_phys_time_ref){
      if(phys_time_ref_frame->registered()){
	phys_time0=phys_time_ref_frame->getRef()[0]-J2000day0;
	cout<<"Trajectory::setup: Using JDframe yields phys_time0="<<phys_time0<<" since start of J2000 epoch. "<<endl;
      }
      else{ //working without a defined frame in this case
	cout<<"Trajectory::setup: Defining JDframe using phys_time0="<<phys_time0<<" since start of J2000 epoch. "<<endl;
	vector<double> ref(1);
	ref[0]=phys_time0+J2000day0;
	phys_time_ref_frame->setRegister(ref);
      }
    }
    *optValue("ref_time")>>tpass_ref_time;
    if(optSet("remap_r0")){
      *optValue("remap_r0_ref_val")>>r0_ref;
      do_remap_r0=true;
    }
    double max_r0;
    *optValue("max_r0")>>max_r0;
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
    
    double r0s=max_r0;
    if(allow_neg_r0)r0s*=2.0;//restricting to positive only (not appropriate with parallax)
    if(do_remap_r0){
      if(allow_neg_r0)cout<<"Trajectory::setup:  Warning remap_r0 not appropriate with ParallaxTrajectory!"<<endl;
      r0s=1.0;
    }
    const int uni=mixed_dist_product::uniform, gauss=mixed_dist_product::gaussian, pol=mixed_dist_product::polar; 
    valarray<double>    centers((initializer_list<double>){ r0s/2.0, tE_max/2,  t0     });
    valarray<double> halfwidths((initializer_list<double>){ r0s/2.0, tE_max/2,  twidth });
    valarray<int>         types((initializer_list<int>){        uni,      uni,  gauss  });
    if(allow_neg_r0 and not do_remap_r0)centers[0]=0;
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
  ///[To be fixed.  As Jeremy points out impact parameter is essentially a 1-D draw, not 2-D, and should thus scale linearly not with cross-section]
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
    addOption("max_r0","Max r0 impact parameter range. (linear prior only,set negative to force positve r0, default=6.0)","6.0");
    addOption("remap_r0","Use remapped r0 coordinate.");
    addOption("remap_r0_ref_val","Cutoff scale for remap of r0.","2.0");
    addOption("tE_max","Uniform prior max in tE. Default=150.0/","150.0");
    addOption("ref_time","Interpret tpass parameter as offset from this time. (Default 0.0","0.0");
  };

};
/////////////////////////////////////////////////////////////////////////////////////////////////////


///This is a ParallaxTrajectory class for when the observer motion relative to the Sun cannot be neglected.
///The nominal version of this class includes the Earth's trajectory accurate to
///1e4 km or about a few minutes motion.  The class is designed to be easily generalizable
///to other orbits by overloading the function that provides the orbital trajectory.
///
///Parameter comments:
/// The SSB velocity is assumed to be constant, and oriented in the "x-direction"
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
    allow_neg_r0=true;
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
  ///First argument should be frame-time
  virtual void get_obs_pos_ssb(double t, double &x, double &y, double &z)const;
  ///Define the observer trajectory velocity in SSB coordinates
  ///For this class we define an approx Earth-Moon barycenter trajectory but that can be overloaded.
  virtual void get_obs_vel_ssb(double t, double &x, double &y, double &z)const;
  ///Using the approx location of the source, transform from ssb in source-lens line-of-sight frame.
  virtual Point ssb_los_transform(double x, double y, double z)const;
  ///Compute observer position offset from ssb in source-lens line-of-sight frame.
  ///Argument takes frame-time
  virtual Point get_obs_pos_offset(double t)const{
    double x=0,y=0,z=0;
    get_obs_pos_ssb(t,x,y,z);
    Point result=ssb_los_transform(x,y,z)*piE;
    if(verbose)
#pragma omp critical 
      {
	cout<<"ParallaxTrajectory::get_obs_pos_offset: x,y,z= "<<x<<","<<y<<","<<z<<endl;
	cout<<"ParallaxTrajectory::get_obs_pos_offset: t= "<<t<<endl;
	cout<<"ParallaxTrajectory::get_obs_pos_offset/norm: = ("<<result.x/piE<<","<<result.y/piE<<")"<<endl;
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
  void equatorial2ecliptic(double source_ra, double source_dec, double &source_lat, double &source_lon){
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
    if(0)
#pragma omp critical 
      {
	cout<<"ParallaxTrajectory::ssb_los_transform: (cos,sin) of source lon="<<source_lon<<" --> ("<<cos_source_lon<<","<<sin_source_lon<<")"<<endl;
	cout<<"ParallaxTrajectory::ssb_los_transform: (cos,sin) of source lat="<<source_lat<<" --> ("<<cos_source_lat<<","<<sin_source_lat<<")"<<endl;
      }
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
