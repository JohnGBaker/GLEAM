//Gravitational lens equation for microlensing
//Written by John G Baker NASA-GSFC (2014-2016)

//#include "glens.hh"
#include "trajectory.hh"
using namespace std;

bool Trajectory::verbose=false;

//
// ******************************************************************
// Point routines **********************************************
// ******************************************************************
//

Point operator+(const Point &a, const Point &b){return Point(a.x+b.x,a.y+b.y);};
Point operator-(const Point &a, const Point &b){return Point(a.x-b.x,a.y-b.y);};
Point operator*(const Point &a, const double &val){return Point(a.x*val,a.y*val);};

//
// ******************************************************************
// Trajectory routines **********************************************
// ******************************************************************
//

// ParallaxTrajectory ***********************************************


///Define the observer trajectory position in SSB coordinates.
///
///We rely on the expression for the Earth-Moon barycenter (EMB) position from
/// http://ssd.jpl.nasa.gov/?planet_pos which are accurate to within a few
/// thousand km.  Note the distance from a point on the Earth or LEO to EMB is 
/// within +/-10000km, a similar scale.  With Earth's orbit at ~30km/s this 
/// level of precision should be compatible with timing precision to within a 
/// few minutes.  For timing, the 1 AU light crossing time is also a few 
/// minutes, so we can similarly ignore that.  If precision better than a few 
/// minutes is relevant then we will need to take these factors into account.
/// The results are returned in the argument in Cartesian SSB coords consistent
/// with ecliptic sky coordinates.
///
void ParallaxTrajectory::get_obs_pos_ssb(double t, double & x, double &y, double &z)const{
  t=t+phys_time0;
  ///We take t to be terrestrial time, TT, time in days since the beginning of 
  ///the J2000 epoch, meaning seconds since J2000 (ts) divided by 86400. 
  /// I.e. t=ts/86400.  Note that the J2000 epoch reference corresponds to 
  /// 2000 Jan 1 11:58:55.816 UTC.  The JPL ephemerides are referenced to 
  ///T_{eph} which is effectively the same as the IAU's TDB.  TDB is nearly 
  ///the same as IAU's terrestrial time TT, corrected for periodic relativistic
  ///changes in the Earth's surface time as the Earth moves in orbit.  The 
  ///corrections are less than 2 ms, so we approximate T_{eph} as TT here. 

  ///Our position comes from JPL's approximate current epoch orbital elements 
  ///( semi-maj axis a, eccen. e, inclination I, mean long. L, long. at
  /// perihelion lonp and long. of ascending node (=0 for earth)) and their 
  /// per-century rates of change. These are:
  const double emax=0.99;//anyway this would be an extremely non-physical time...
  const double degrad=180.0/M_PI;
  const double a0=1.00000261,e0=0.01671123,I0=-0.00001531/degrad,
    L0=100.46457166/degrad,lonp0=102.93768193/degrad;
  const double adot=5.62e-6,edot=-43.92e-6,Idot=-0.01294668/degrad,
    Ldot=35999.37244981/degrad,lonpdot=0.32327364/degrad;
  //inst. values
  double Tcen=t/36525.;
  double a=a0+adot*Tcen,e=e0+edot*Tcen,I=I0+Idot*Tcen,L=L0+Ldot*Tcen,lonp=lonp0+lonpdot*Tcen;
  if(e<-emax)e=-emax;
  if(e>emax)e=emax;
  //argument of perihelion argp=lonp here
  //mean anomaly M (in range -pi<=M<PI, eccentric anom. E (solves M=E-e*sin(E),
  double M=(floor(((L-lonp)/M_PI+1)/2.0)-1.0)*M_PI*2,E=M+e*sin(M),Eold=0;
  //solve for E
  while(abs(E-Eold)>1e-13*abs(Eold)){Eold=E;E=Eold+(M-Eold+e*sin(Eold))/(1-e*cos(Eold));}
  //perihl. cartesian coords
  double xper=a*(cos(E)-e),yper=a*sqrt(1-e*e)*sin(E);
  double cp=cos(lonp),sp=sin(lonp),cI=cos(I),sI=sin(I);
  if(verbose)
#pragma omp critical 
    {
      cout<<"ParallaxTrajectory::get_obs_pos_ssb: cp,xper,sp,yper,E,Eold:"<<cp<<", "<<xper<<", "<<sp<<", "<<yper<<", "<<E<<", "<<Eold<<endl;
      cout<<"L,e,M,a:"<<L<<", "<<e<<", "<<M<<", "<<a<<endl;
      cout<<"phys_time0="<<phys_time0<<endl;
    }
  //ecliptic coords
  x=cp*xper-sp*yper;
  double rhoyz=(sp*xper+cp*yper);
  y=rhoyz*cI;
  z=rhoyz*sI;
};

///Here we just compute the second order numerical derivative.  That is probably good enough...but we could do analytic..
///It is not hard to write down an analytic version of this to be computed along-with the position computation. Basically
///you can certainly treat a and e as constant, and you can probably also treat the other elements except for M,E as constant
///in taking the derivatives.  The trouble is you have to do this computation more or less along with the position computation.
///We could just copy the above computation here with the changes and recompute those derivs.  We could also compute and
///save velocities along with the position computation, and then chack and don't recompute if the evaluation time is identical.
///That would save up to a factor of almost 3 over this computation if velocity and position are always evals consequtively at
///the same points, or a factor of 2 if not.
void ParallaxTrajectory::get_obs_vel_ssb(double t, double & rdot, double &thdot, double &phdot)const{
  double rp,thp,php,rm,thm,phm;
  double delta=1.0; //1 day is small
  get_obs_pos_ssb(t+delta,rp,thp,php);
  get_obs_pos_ssb(t-delta,rm,thm,phm);
  rdot=(rp-rm)/(2*delta);
  thdot=(thp-thm)/(2*delta);
  phdot=(php-phm)/(2*delta);
}; 
  ///Using the approx location of the source, transform from ssb in source-lens line-of-sight frame.
Point ParallaxTrajectory::ssb_los_transform(double x, double y, double z)const{
  ///Here we must assume some coordinate frame orientation for the observer plane. 
  ///In particular we assume that the ecliptic intersects the observer plane along the x-axis
  ///and SSB north projects onto the y-axis.  This transformation breaks down if the source-lens
  ///lies on/near one of the SSB poles, as is roughly the case for the LMC.  An alternative would
  /// The transformation from BBS (r,\theta,\phi) to the observation plane x-y-z coordinates is then:
  ///   \f$ x = - r \sin\theta \sin(\phi-l) \f$
  ///   \f$ y = - r \sin\theta \cos(\phi-l) \sin b + r \cos\theta \cos b \f$
  ///   \f$ z = - r \sin\theta \cos(\phi-l) \cos b - r \cos\theta \sin b \f$
  ///with l=source_lon and b=source_lat.  We have included z for completeness, though it is irrelevant.
  ///In terms of SSB cartesian coordinates: \f$ z_{SSB} = r\cos\theta, x_{SSB} = r\sin\theta\cos\phi, etc. \f$ 
  ///So we get:
  ///   xnew = -y*cl + x*sl
  ///   ynew = -(x*cl+y*cl)*sb +z*cb
  ///   znew = -(x*cl+y*cl)*cb -z*sb
  /// where cl=cos(source_lon), etc.
  double xnew,ynew;
  xnew = -y*cos_source_lon + x*sin_source_lon;
  ynew = -(x*cos_source_lon + y*sin_source_lon)*sin_source_lat + z*cos_source_lat;
  //last we rotate to lens frame
  double xlens,ylens;
  xlens = xnew * cos_dphi - ynew * sin_dphi;
  ylens = xnew * sin_dphi + ynew * cos_dphi;
  if(verbose)
#pragma omp critical 
    {
      cout<<"ParallaxTrajectory::ssb_los_transform: xnew,ynew= "<<xnew<<","<<ynew<<endl;
      cout<<"ParallaxTrajectory::ssb_los_transform: xlens,ylens= "<<xlens<<","<<ylens<<endl;
    }
  
  return Point(xlens,ylens);
};


