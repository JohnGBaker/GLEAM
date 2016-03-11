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

#endif
