//Gravitational lens equation for microlensing
//Written by John G Baker NASA-GSFC (2014)

#ifndef GLENS_HH
#define GLENS_HH
#include <vector>
#include <iostream>
#include <cmath>
#include <sstream>
#include <iomanip>
#include "bayesian.hh"
using namespace std;
extern bool debug;

///Code for computing binary microlensing light curves and images.
///This part deals with the binary lens and the lens map.
///Generally the reference frame is defined by the source and lens centers
///with the observer moving through the some observing plane.

///We are thinking here of a binary lens, but there are other possibilities.
///The core methods for any lens are the lensing map, magnification, etc.
///Finite source effects would generally require integration over the source 
///plane.
///There are various possible implementations for a point lens everything can 
///be done analytically.  For a fixed binary, the lens map can be reduced to a 
///polynomial root-finding problem.  Generically, we can grid the lens plane and
///apply brute force. This would need to be done only once per set of lens params (masses,separation).

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
class Trajectory {
  //Transition notes:
  //We are getting ready to promote Trajectory to a bayes_component, and then
  //to move definition of trajectory-related parameters into here.
  //My current thinking is that everything is done nominally in a frame with
  //the source-CM to lens-CM line-of-sight as origin in the observer plane.
  //Trajectory parameters indicate where the Sun/Earth passes through that
  //plane, thus Trajectory parameters include at least r0, [phi?], and probably
  //also tmax.  An issue is that, as of now, tmax is defined in physical
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
  double tE,tmax;
  double r0,phi;
public:
  Trajectory(Point pos0,Point vel0,double t_end,double cadence, double ts=0):p0(pos0),v0(vel0),cad(cadence),tf(t_end),ts(ts),toff(0),have_times(false),tE(1),tmax(0){
    //cout<<"Trajectory({"<<pos0.x<<","<<pos0.y<<"},{"<<vel0.x<<","<<vel0.y<<"},...)"<<endl;
  };
  Trajectory(Point pos0,Point vel0):p0(pos0),v0(vel0),cad(1),tf(1),ts(0),have_times(false),tE(1),tmax(0){
    //cout<<"Trajectory({"<<pos0.x<<","<<pos0.y<<"},{"<<vel0.x<<","<<vel0.y<<"})"<<endl;
  };
  Trajectory():Trajectory(Point(0,0),Point(1,0)){};
  virtual ~Trajectory(){};//Need virtual destructor to allow derived class objects to be deleted from pointer to base.
  virtual Trajectory* clone(){
    return new Trajectory(*this);
  };
  virtual void setup(const Point &pos0, const Point &vel0){p0=pos0;v0=vel0;};
  ///Set the required eval times.  (Times are "phys" times, but "phys"=frame if tE=1,tmax=0)
  virtual void set_times(vector<double> times,double toff){this->times=times;ts=times[0];tf=times.back();have_times=true;this->toff=toff;};
  ///Return start time. (Times are "phys" times, but "phys"=frame if tE=1,tmax=0)
  virtual double t_start()const {return ts;};
  ///Return end time. (Times are "phys" times, but "phys"=frame if tE=1,tmax=0)
  virtual double t_end()const {return tf;};
  virtual int Nsamples()const {if(have_times)return times.size(); else return (int)((t_end()-t_start())/cad)+1;};
  virtual double get_obs_time(int ith)const {if(have_times)return times[ith]; else return ts-toff+cad*ith;};
  virtual double get_phys_time(double frame_time){return frame_time*tE+tmax;}
  virtual double get_frame_time(double phys_time){return (phys_time-tmax)/tE;}
  ///Argument takes frame time below
  virtual Point get_obs_pos(double t)const {double x=p0.x+(t-toff)*v0.x,y=p0.y+(t-toff)*v0.y;return Point(x,y);};
  virtual Point get_obs_vel(double t)const {return v0;};
  virtual string print_info()const {ostringstream s;s<<"Trajectory({"<<p0.x<<","<<p0.y<<"},{"<<v0.x<<","<<v0.y<<"})"<<endl;return s.str();};
};

///Here is ParallaxTrajectory class for when the observer motion cannot be neglected.
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
/// So far, the trajectory parameters are interpreted in the lens-frame, with an angle parameter
/// which specifies the rotation between the lens frame and the velocity direction of the
/// Earth-sun projected to the image plane.  Since the velocity is assumed to be constant, we
/// are free to more precisely interpret this angle as orientation at the point of closest
/// approach of the Sun's path to the COM line-of-sight origin.
/// To complete the specification we need another angle relating the Sun's path to the physical
/// frame (e.g. of the ecliptic). That additional angle is a heliocentric version of the so-called
/// parallax direction (See Gould2014[arxiv:1408.0797] for somewhat relevant discussion). To give
/// this parameter a name, based on that literature, we call it \f$ \phi_\mu \f$ (ie phimu). We also
/// need a parameter relating the scale of 1AU, the lenght unit for our ephemeris, to the Einstein
/// units in the image plane.  By convention (see eg Gaudi2010[arxiv[1002:0332]) this parameter is
/// called \f$ \pi_E \f$. A 1AU displacement in the observer plane corresponds to a piE displacement
/// in \beta in Einstein units.  This is related to the mass and distances by:
/// piE^2 = (1AU/Drel/thetaE)^2 = (1/4)(1AU/Drel)(1AU/M)
/// where Drel = D_s*D_l/(D_s-D_l)
class ParallaxTrajectory : public Trajectory {
protected: 
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
  ParallaxTrajectory(Point pos0, Point vel0, double t_end, double caden, double source_ra, double source_dec, double t2000off, double ts=0):Trajectory(pos0,vel0,t_end,caden,ts){
    equatorial2ecliptic(source_ra,source_dec,source_lat,source_lon);
    cos_source_lon=cos(source_lon);
    sin_source_lon=sin(source_lon);
    cos_source_lat=cos(source_lat);
    sin_source_lat=sin(source_lat);
  };
  virtual ParallaxTrajectory* clone(){
    return new ParallaxTrajectory(*this);
  };
  virtual void setup(double pieE,double phimu){};//2TRAJ: put something here
  virtual void rotate(double pieE,double phimu){};//2TRAJ: put something here
  Point get_obs_pos(double t)const {return Trajectory::get_obs_pos(t)+get_obs_pos_offset(t);};
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
    double x,y,z;
    get_obs_pos_ssb(t,x,y,z);
    return ssb_los_transform(x,y,z)*piE;
  };
  ///Compute observer velocity offset from ssb in source-lens line-of-sight frame.
  virtual Point get_obs_vel_offset(double t)const{
    double x,y,z;
    get_obs_vel_ssb(t,x,y,z);
    return ssb_los_transform(x,y,z)*piE;
  };
  ///approximately convert equatorial to ecliptic coordinates
  void equatorial2ecliptic(double source_ra, double source_dec, double source_lat, double source_lon){
    const double eps=23.4372; //as of 2016.0, decreasing at 0.00013/yr
    const double ceps=cos(eps),seps=sin(eps);
    double sa=sin(source_ra),ca=cos(source_ra);
    double sd=sin(source_dec),cd=cos(source_dec);
    source_lon=atan2(sa*ceps+sd/cd*seps,ca);
    source_lat=sd*ceps-cd*seps*sa;
  };
};

///This is a generic abstract base class for thin gravitational lens objects.
class GLens :public bayes_component{
protected:
  int NimageMax;
  static const double constexpr dThTol=1e-9;
  ///For temporary association with a trajectory;
  const Trajectory *trajectory;
  ///Transform from trajectory frame to lens frame
  virtual Point traj2lens(const Point tp)const {return tp;};
  virtual Point lens2traj(const Point tp)const {return tp;};
  virtual Point traj2lensdot(const Point tp)const {return tp;};
public:
  virtual Point get_obs_pos(const Trajectory & traj,double time)const{return traj2lens(traj.get_obs_pos(time));};
  virtual Point get_obs_vel(const Trajectory & traj,double time)const{return traj2lensdot(traj.get_obs_vel(time));};
protected:
  ///static access for use with non-member GSL integration routines
  static int GSL_integration_func (double t, const double theta[], double thetadot[], void *instance);
  static int GSL_integration_func_vec (double t, const double theta[], double thetadot[], void *instance);
  //static Point get_obs_pos(const GLens* instance, const Trajectory & traj,double time){return instance->get_obs_pos(traj,time);};
  //static Point get_obs_vel(const GLens* instance, const Trajectory & traj,double time){return instance->get_obs_vel(traj,time);};
  ///For use with GSL integration
  double kappa=.1;
  int Ntheta;
  bool use_integrate,have_integrate,do_verbose_write;
  double GL_int_tol,GL_int_mag_limit;
  virtual bool testWide(const Point & p,double scale)const{return false;};//test conditions to revert to perturbative inversion
public:
  virtual ~GLens(){};//Need virtual destructor to allow derived class objects to be deleted from pointer to base.
  GLens(){have_integrate=false;do_verbose_write=false;};
  virtual GLens* clone(){return new GLens();};
  ///Lens map: map returns a point in the observer plane from a point in the lens plane.
  virtual Point map(const Point &p){cout<<"GLens::map: This should be a single lens of unit mass. It's a simple function: place it here if you need it."<<endl;exit(1);};
  ///Inverse sens map: invmap returns a set of points in the lens plane which map to some point in the observer plane.  Generally multivalued;
  virtual vector<Point> invmap(const Point &p){cout<<"GLens::invmap: This should be a single lens of unit mass. It's a simple function: place it here if you need it."<<endl;exit(1);};
  ///Given a point in the lens plane, return the magnitude
  virtual double mag(const Point &p){cout<<"GLens::mag: This should be a single lens of unit mass. It's a simple function: place it here if you need it."<<endl;exit(1);};
  ///Given a set of points in the lens plane, return the combined magnitude
  virtual double mag(const vector<Point> &plist){
    double m=0;
    for(Point p : plist){
      m+=abs(mag(p));//syntax is C++11
      if(debug)cout<<"    ("<<p.x<<","<<p.y<<") --> mg="<<m<<endl;
    }
    if(plist.size()==0)return 1; // to more gracefully fail in trivial regions
    return m;
  };
  ///returns J=det(d(map(p))/dp)^-1, sets, j_ik = d(map(pi))/dpk
  virtual double jac(const Point &p,double &j00,double &j01,double &j10,double &j11){cout<<"GLens::jac: This should be a single lens of unit mass. It's a simple function: place it here if you need it."<<endl;exit(1);};
  ///returns J=det(d(map(p))/dp))^-1, sets, j_ik = (d(map(pi))/dpk)^-1
  virtual double invjac(const Point &p,double &j00,double &j01,double &j10,double &j11){cout<<"GLens::invjac: This should be a single lens of unit mass. It's a simple function: place it here if you need it."<<endl;exit(1);};;
  ///compute images and magnitudes along some trajectory
  void compute_trajectory (const Trajectory &traj, vector<double> &time_series, vector<vector<Point> > &thetas_series, vector<int> &index_series,vector<double>&mag_series,bool integrate=false);
  void set_integrate(bool integrate_or_not){use_integrate=integrate_or_not;have_integrate=true;}
  //For the Optioned interface:
  virtual void addOptions(Options &opt,const string &prefix="");
  void static addStaticOptions(Options &opt){
    GLens l;
    l.addTypeOptions(opt);
  };
  void addTypeOptions(Options &opt){
    Optioned::addOptions(opt,"");
    addOption("binary_lens","Apply a binary lens model.");
  };  
  virtual void setup();
  virtual string print_info()const{ostringstream s;s<<"GLens()"<<(have_integrate?(string("\nintegrate=")+(use_integrate?"true":"false")):"")<<endl;return s.str();};
  //For stateSpaceInterface
  virtual void defWorkingStateSpace(const stateSpace &sp){haveWorkingStateSpace();};
  virtual stateSpace getObjectStateSpace()const{return stateSpace(0);};
  virtual void setState(const state &s){checkWorkingStateSpace();};
;
  //getCenter provides *trajectory frame* coordinates for the center. Except for with -2, which give the lens frame CM. 
  virtual Point getCenter(int option=-2)const{return Point(0,0);};
  //Write a magnitude map to file.  
  //Points in this function and its arguments are in *trajectory frame* coordinates 
  virtual void writeMagMap(ostream &out, const Point &LLcorner,const Point &URcorner,int samples){//,bool output_nimg=false){
    cout<<"GLens::writeMagMap from ("<<LLcorner.x<<","<<LLcorner.y<<") to ("<<URcorner.x<<","<<URcorner.y<<")"<<endl;
    double dx=(URcorner.x-LLcorner.x)/(samples-1);    
    double dy=(URcorner.y-LLcorner.y)/(samples-1);    
    //cout<<"mag-map ranges from: ("<<x0<<","<<y0<<") to ("<<x0+width<<","<<y0+width<<") stepping by: "<<dx<<endl;
    int output_precision=out.precision();
    ios_base::fmtflags flags=out.flags();
    //cout<<"writeMagMap:output_precision="<<output_precision<<endl;
    double ten2prec=pow(10,output_precision-2);
    out<<"#x  y  magnification"<<endl;
    for(double y=LLcorner.y;y<=URcorner.y;y+=dy){
      //Trajectory traj(Point(LLcorner.x+xc,y), Point(1,0), URcorner.x-LLcorner.x, dx);
      Trajectory traj(Point(LLcorner.x,y), Point(1,0), URcorner.x-LLcorner.x, dx);
      vector<int> indices;
      vector<double> times,mags;
      vector<vector<Point> >thetas;
      compute_trajectory(traj,times,thetas,indices,mags);//2TRAJ:Check that the new interface doesn't break this.
      for(int i : indices){
	//Point b=traj.get_obs_pos(times[i]);//TRAJ:Allow non-trivial transformation from observer-plane coords frame to lens frame coords
	Point b=traj.get_obs_pos(times[i]);//we want the result in traj frame, to match the dump_trajectory output
	double mtruc=floor(mags[i]*ten2prec)/ten2prec;
	//out.precision(output_precision);
	out<<b.x<<" "<<b.y<<" "<<setiosflags(ios::scientific)<<mtruc<<resetiosflags(flags);
	if(do_verbose_write){
	  out<<" "<<thetas[i].size();
	  if(true){
	    for(int j=0;j<thetas[i].size();j++){
	      out<<" "<<thetas[i][j].x<<" "<<thetas[i][j].y;
	    }	
	  }
	}
	out<<endl;
	
      }
      out<<endl;
    }	  
  };   
  void verboseWrite(bool state=true){do_verbose_write=state;};
  
};

///A rigid binary lens implementation
///Working in units of total mass Einstein radius, the only parameters are
/// mass ratio q and separation L
class GLensBinary : public GLens{
  double q;
  double L;
  double phi0,sin_phi0,cos_phi0;
  Point cm;//center of mass in lens frame
  //mass fractions
  double nu;
  vector<Point> invmapAsaka(const Point &p);
  vector<Point> invmapWittMao(const Point &p);
  double rWide;
  //parameter handling
  double q_ref;
  bool do_remap_q;
  int idx_q,idx_L,idx_phi0;
  virtual Point traj2lens(const Point tp)const {
    //Note: this looks like a -phi0 rotation here because phi0 is the rotation from the observer to lens *frame axis* and here we transform the ccordinates
    return Point(cm.x+tp.x*cos_phi0-tp.y*sin_phi0,cm.y+tp.x*sin_phi0+tp.y*cos_phi0);};
  virtual Point lens2traj(const Point tp)const {
    return Point((tp.x-cm.x)*cos_phi0+(tp.y-cm.y)*sin_phi0,-(tp.x-cm.x)*sin_phi0+(tp.y-cm.y)*cos_phi0);};
  virtual Point traj2lensdot(const Point tp)const {
    return Point(tp.x*cos_phi0-tp.y*sin_phi0,cm.y+tp.x*sin_phi0+tp.y*cos_phi0);};
  ///test conditions to revert to perturbative inversion
  bool testWide(const Point & p,double scale)const{
    double rs=rWide*scale;
    if(rs<=0)return false;
    double r2=p.x*p.x+p.y*p.y;
    //if(( L>rs||r2>rs*rs)&&scale!=1.0)cout<<sqrt(r2)<<" <> "<<rs<<" <> "<<L<<endl;
    //return L>rs||r2>rs*rs;
    return L>rs||r2>rs*rs||(q+1/q)>2*rs*rs;
  };  
public:
  GLensBinary(double q=1,double L=1,double phi0=0);
  virtual GLensBinary* clone(){
    return new GLensBinary(*this);
  };
  virtual void setup();
  Point map(const Point &p);
  //For the GLens interface:
  vector<Point> invmap(const Point &p);
  vector<Point> invmapWideBinary(const Point &p);
  double mag(const Point &p);
  using  GLens::mag;
  double jac(const Point &p,double &j00,double &j01,double &j10,double &j11);
  double invjac(const Point &p,double &j00,double &j01,double &j10,double &j11);
  //specific to this class:
  double get_q(){return q;};
  double get_L(){return L;};
  double set_WideBinaryR(double r){rWide=r;};
  virtual string print_info()const{ostringstream s;s<<"GLensBinary(q="<<q<<",L="<<L<<")"<<(have_integrate?(string("\nintegrate=")+(use_integrate?"true":"false")):"")<<endl;return s.str();};
  //virtual string print_info()const{ostringstream s;s<<"GLensBinary(q="<<q<<",L="<<L<<")"<<endl;return s.str();};

  ///From StateSpaceInterface (via bayes_component)
  ///
  void defWorkingStateSpace(const stateSpace &sp){
    checkSetup();
    if(do_remap_q)idx_q=sp.requireIndex("s(1+q)");
    else idx_q=sp.requireIndex("logq");
    idx_L=sp.requireIndex("logL");
    idx_phi0=sp.requireIndex("phi0");
    haveWorkingStateSpace();
  };  
  ///Set up the output stateSpace for this object
  ///
  stateSpace getObjectStateSpace()const{
    checkSetup();//Call this assert whenever we need options to have been processed.
    stateSpace space(3);
    string names[]={"logq","logL","phi0"};
    if(do_remap_q)names[3]="s(1+q)";
    space.set_bound(2,boundary(boundary::wrap,boundary::wrap,0,2*M_PI));//set 2-pi-wrapped space for phi0.
    space.set_names(names);  
    return space;
  };
  virtual void addOptions(Options &opt,const string &prefix=""){
    GLens::addOptions(opt,prefix);
    addOption("remap_q","Use remapped mass-ratio coordinate.");//2TRAJLENS:move to GLensBinary
    addOption("q0","Prior max in q (with q>1) with remapped q0. Default=1e5/","1e5");//2TRAJLENS:move to GLensBinary
  };
  ///Set state parameters
  ///
  void setState(const state &st){;
    //  logL separation (in log10 Einstein units)
    //  q mass ratio
    //  alignment angle phi0, (binary axis rel to trajectory direction) at closest approach point
    checkWorkingStateSpace();
    double f_of_q=st.get_param(idx_q);//either log_q or remapped q
    double logL=st.get_param(idx_L);
    phi0=st.get_param(idx_phi0);

    L=pow(10.0,logL);
    //(See discussion in remap_q() above.) otherwise sofq==ln_q.
    if(do_remap_q)q=-1.0+(q_ref+1.0)/sqrt(1.0/f_of_q-1.0);
    else q=pow(10.0,f_of_q);
    cos_phi0=cos(phi0);
    sin_phi0=sin(phi0);
    nu=1/(1+q);
    cm=Point((q/(1.0+q)-0.5)*L,0);
  };
  ///This is a class-specific variant which (probably no longer needed)
  void setState(double q_, double L_){;
    checkWorkingStateSpace();
    q=q_;
    L=L_;
    nu=1/(1+q);
    cm=Point((q/(1.0+q)-0.5)*L,0);
  };

  //getCenter provides *trajectory frame* coordinates for the center.  //be consistent with traj2lens
  Point getCenter(int option=-2)const{
    double x0=0,y0=0,width=0,wx,wy;
    //cout<<"q,L,option="<<q<<", "<<L<<", "<<option<<endl;
    //center on {rminus-CoM,CoM-CoM,rplus-CoM}, when cent={-1,0,1} otherwise CoM-nominalorigin;
    switch(option){
    case -1://minus lens rel to CoM
      x0=-0.5*L;
      break;
    case 0:
      x0=0;
      break;
    case 1:
      x0=0.5*L;//plus lens rel to CoM
      break;
    default:
      x0=cm.x;
      //cout<<"case def"<<endl;
    }
    //cout<<" GLensBinary::getCenter("<<option<<"):Returning x0="<<x0<<endl;
    //This returns the result in lens frame
    return lens2traj(Point(x0,0));
  };
  
protected:
  
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

};

#endif
