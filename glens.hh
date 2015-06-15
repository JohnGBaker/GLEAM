//Gravitational lens equation for microlensing
//Written by John G Baker NASA-GSFC (2014)

#ifndef GLENS_HH

#include <vector>
#include <iostream>
#include <cmath>
#include <sstream>
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
  Point(double x, double y):x(x),y(y){};
}Point;
//Point minus(const Point &a, const Point &b){return Point(a.x-b.x,a.y-b.y)};

///Next is a class for trajectories through the observer plane
///base class implements a straight-line trajectory
class Trajectory {
protected:
  Point p0;
  Point v0;
  double t0;
  double tf;
  double cad;
  double toff;
  bool have_times;
  vector<double> times;
public:
  Trajectory(Point pos0,Point vel0,double t_end,double cadence, double t0=0):p0(pos0),v0(vel0),cad(cadence),tf(t_end),t0(t0),toff(0),have_times(false){
    //cout<<"Trajectory({"<<pos0.x<<","<<pos0.y<<"},{"<<vel0.x<<","<<vel0.y<<"},...)"<<endl;
  };
  Trajectory(Point pos0,Point vel0):p0(pos0),v0(vel0),cad(1),tf(1),t0(0),have_times(false){
    //cout<<"Trajectory({"<<pos0.x<<","<<pos0.y<<"},{"<<vel0.x<<","<<vel0.y<<"})"<<endl;
  };
  ///set the required eval times.
  virtual void set_times(vector<double> times,double toff){this->times=times;t0=times[0];tf=times.back();have_times=true;this->toff=toff;};//cout<<"tf="<<tf<<", times="<<this->times.size()<<", times["<<times.size()-1<<"]="<<times[times.size()-1]<<endl;};
  virtual double t_start()const {return t0;};
  virtual double t_end()const {return tf;};
  virtual double get_obs_time(int ith)const {if(have_times)return times[ith]; else return t0-toff+cad*ith;};
  virtual Point get_obs_pos(double t)const {double x=p0.x+(t-toff)*v0.x,y=p0.y+(t-toff)*v0.y;return Point(x,y);};
  virtual Point get_obs_vel(double t)const {return v0;};
  virtual string print_info()const {ostringstream s;s<<"Trajectory({"<<p0.x<<","<<p0.y<<"},{"<<v0.x<<","<<v0.y<<"})"<<endl;return s.str();};
};

///Next is a class for trajectories through the observer plane
///class implements a straight-line trajectory
class ParallaxTrajectory : public Trajectory {
public:
  ParallaxTrajectory(Point pos0,Point vel0, double t_end, double caden, double au, double lat,double lon, double psi0, double t0=0):Trajectory(pos0,vel0,t_end,caden,t0){};
  double t_start()const {return t0;};
  double t_end()const {return tf;};
  double get_obs_time(int ith)const {return t0+cad*ith;};
  Point get_obs_pos(double t)const {double x=p0.x+t*v0.x,y=p0.y+t*v0.y;return Point(x,y);};
  Point get_obs_vel(double t)const {return v0;};
};

///This is a generic abstract base class for thin gravitational lens objects.
class GLens {
protected:
  int NimageMax;
  static const double constexpr dThTol=1e-9;
  ///For temporary association with a trajectory;
  const Trajectory *trajectory;
  ///For use with GSL integration
  static int GSL_integration_func (double t, const double theta[], double thetadot[], void *instance);
  static int GSL_integration_func_vec (double t, const double theta[], double thetadot[], void *instance);
  ///For use with GSL integration
  static const int kappa=.1;
  int Ntheta;
public:
  ///Lens map: map returns a point in the observer plane from a point in the lens plane.
  virtual Point map(const Point &p)=0;
  ///Inverse sens map: invmap returns a set of points in the lens plane which map to some point in the observer plane.  Generally multivalued;
  virtual vector<Point> invmap(const Point &p)=0;
  ///Given a point in the lens plane, return the magnitude
  virtual double mag(const Point &p)=0;
  ///Given a set of points in the lens plane, return the combined magnitude
  virtual double mag(const vector<Point> &plist){
    double m=0;
    for(Point p : plist){
      m+=abs(mag(p));//syntax is C++11
      if(debug)cout<<"    ("<<p.x<<","<<p.y<<") --> mg="<<m<<endl;
    }
    if(plist.size()==0)return 1; //hack to more gracefully fail in trivial regions
    return m;
  };
  ///returns J=det(d(map(p))/dp)^-1, sets, j_ik = d(map(pi))/dpk
  virtual double jac(const Point &p,double &j00,double &j01,double &j10,double &j11)=0;
  ///returns J=det(d(map(p))/dp))^-1, sets, j_ik = (d(map(pi))/dpk)^-1
  virtual double invjac(const Point &p,double &j00,double &j01,double &j10,double &j11)=0;
  ///compute images and magnitudes along some trajectory
  void compute_trajectory (const Trajectory &traj, vector<double> &time_series, vector<vector<Point> > &thetas_series, vector<int> &index_series,vector<double>&mag_series,bool integrate=false);
};

///A rigid binary lens implementation
///Working in units of total mass Einstein radius, the only parameters are
/// mass ratio q and separation L
class GLensBinary : public GLens{
  double q;
  double L;
  //mass fractions
  double nu;
  vector<Point> invmapAsaka(const Point &p);
  vector<Point> invmapWittMao(const Point &p);
  double rWide;
public:
  GLensBinary(double q,double L);
  Point map(const Point &p);
  vector<Point> invmap(const Point &p);
  vector<Point> invmapWideBinary(const Point &p);
  double mag(const Point &p);
  using  GLens::mag;
  double jac(const Point &p,double &j00,double &j01,double &j10,double &j11);
  double invjac(const Point &p,double &j00,double &j01,double &j10,double &j11);
  double get_q(){return q;};
  double get_L(){return L;};
  double set_WideBinaryR(double r){rWide=r;};
};

#define GLENS_HH
#endif
