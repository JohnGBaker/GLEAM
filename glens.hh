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
#include "trajectory.hh"
#include <complex>

using namespace std;
extern bool debug;

///Code for computing microlensing light curves and images.
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

///This is a generic (abstract) base class for thin gravitational lens objects.
class GLens :public bayes_component{
protected:
  int NimageMax;
  int NimageMin;
  static const double constexpr dThTol=1e-9;
  ///finite_source
  bool do_finite_source;
  int finite_source_method;
  int idx_log_rho_star;
  double source_radius;
  //StateSpace and Prior
  stateSpace GLSpace;

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
  virtual int poly_root_integration_func_vec (double t, const double theta[], double thetadot[], void *instance){
    cout<<"GLens::poly_root_integration_func:  Not defined!"<<endl;exit(1);};
  //static Point get_obs_pos(const GLens* instance, const Trajectory & traj,double time){return instance->get_obs_pos(traj,time);};
  //static Point get_obs_vel(const GLens* instance, const Trajectory & traj,double time){return instance->get_obs_vel(traj,time);};
  ///For use with GSL integration
  double kappa=.1;
  int Ntheta;
  bool use_integrate,have_integrate,do_verbose_write;
  double GL_int_tol,GL_int_mag_limit;
  virtual bool testWide(const Point & p,double scale)const{return false;};//test conditions to revert to perturbative inversion
  //Utility for allowing incremental update of nearby solutions.
  bool have_saved_soln;
public:
  virtual ~GLens(){};//Need virtual destructor to allow derived class objects to be deleted from pointer to base.
  GLens(){typestring="GLens";option_name="SingleLens";option_info="Single point-mass lens";have_integrate=false;do_verbose_write=false;have_saved_soln=false;NimageMax=2;NimageMin=2;do_finite_source=false;idx_log_rho_star=-1;};
  virtual GLens* clone(){return new GLens(*this);};
  ///Lens map: map returns a point in the observer plane from a point in the lens plane.
  virtual Point map(const Point &p){
    long double x=p.x,y=p.y,rsq=x*x+y*y,c=(1.0L-1.0L/rsq);;
    //cout<<"map: x,y,c"<<x<<", "<<y<<", "<<c<<endl;
    return Point(x*c,y*c);
  };
  ///Inverse sens map: invmap returns a set of points in the lens plane which map to some point in the observer plane.  Generally multivalued;
  virtual vector<Point> invmap(const Point &p){
    long double x=p.x,y=p.y,rsq=x*x+y*y,c0=sqrt(1.0L+4.0L/rsq);
    //cout<<"map: x,y,r2"<<x<<", "<<y<<", "<<rsq<<endl;
    vector<Point> thetas(2);
    double c=(1.0L+c0)/2.0L;
    thetas[0]=Point(x*c,y*c);
    //cout<<"c0="<<c<<endl;
    c=(1.0L-c0)/2.0L;
    //cout<<"c1="<<c<<endl;
    thetas[1]=Point(x*c,y*c);
    return thetas;
  };
  ///Given a point in the lens plane, return the magnitude
  virtual double mag(const Point &p){
    long double x=p.x,y=p.y,rsq=x*x+y*y,r4=rsq*rsq;
    return 1.0L/(1.0L-r4);
  };
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
  ///Compute the Laplacian of the local image magnification explicitly
  virtual double Laplacian_mu(const Point &p)const;
  ///Compute the complex lens shear, and some number of its derivatives  
  virtual vector<complex<double> > compute_shear(const Point &p, int nder)const;
  ///compute images and magnitudes along some trajectory
  virtual void finite_source_compute_trajectory (const Trajectory &traj, vector<double> &time_series, vector<vector<Point> > &thetas_series, vector<double>&mag_series, ostream *out=NULL);
  void compute_trajectory (const Trajectory &traj, vector<double> &time_series, vector<vector<Point> > &thetas_series, vector<int> &index_series,vector<double>&mag_series,bool integrate=false);
  void inv_map_curve(const vector<Point> &curve, vector<vector<Point> > &curves_images, vector<vector<double>> &curve_mags);
  //Note that the centroid is returned in p, and the variance is returned in var
  static double _image_area_mag_dummy_variance;
  void image_area_mag(Point &p, double radius, int & N, double &magnification, double &var=_image_area_mag_dummy_variance, ostream *out=NULL);
  void set_integrate(bool integrate_or_not){use_integrate=integrate_or_not;have_integrate=true;}
  //For the Optioned interface:
  virtual void addOptions(Options &opt,const string &prefix="");
  /*
  void static addStaticOptions(Options &opt){
    GLens l;
    l.addTypeOptions(opt);
  };
  void addTypeOptions(Options &opt){
    Optioned::addOptions(opt,"");
    addOption("binary_lens","Apply a binary lens model.");
    };*/  
  virtual void setup();
  virtual string print_info()const{ostringstream s;s<<"GLens()"<<(have_integrate?(string("\nintegrate=")+(use_integrate?"true":"false")):"")<<endl;return s.str();};
  //For stateSpaceInterface
  virtual void defWorkingStateSpace(const stateSpace &sp){
    if(do_finite_source)idx_log_rho_star=sp.requireIndex("log_rho_star");
    haveWorkingStateSpace();
  };
  virtual void setState(const state &st){
    bayes_component::setState(st);
    cout<<"idx_log_rho_star="<<idx_log_rho_star<<endl;
    if(do_finite_source)source_radius=pow(10.0,st.get_param(idx_log_rho_star));
    cout<<"source_radius="<<source_radius<<endl;
  };
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
      Trajectory traj(Point(LLcorner.x,y), Point(1,0), URcorner.x-LLcorner.x, dx);
      vector<int> indices;
      vector<double> times,mags;
      vector<vector<Point> >thetas;
      compute_trajectory(traj,times,thetas,indices,mags);//2TRAJ:Check that the new interface doesn't break this.
      for(int i : indices){
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
  shared_ptr<const sampleable_probability_function> parentPrior;//remember parent prior so we can replace nativePrior
  shared_ptr<const sampleable_probability_function> binaryPrior;
  stateSpace GLBinarySpace;
  double q;
  double L;
  double phi0,sin_phi0,cos_phi0;
  Point cm;//center of mass in lens frame
  //mass fractions
  double nu;
  vector<Point> invmapAsaka(const Point &p);
  vector<Point> invmapWittMao(const Point &p,bool no_check=false);
  virtual int poly_root_integration_func_vec (double t, const double theta[], double thetadot[], void *instance);
  //complex<double> saved_roots[6];
  vector<Point> theta_save;
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
  ///returns J=det(d(map(p))/dp)^-1, sets, j_ik = d(map(pi))/dpk
  double jac(const Point &p,double &j00,double &j01,double &j10,double &j11);
  double invjac(const Point &p,double &j00,double &j01,double &j10,double &j11);
  vector<complex<double> > compute_shear(const Point &p, int nder)const;
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
    GLens::defWorkingStateSpace(sp);
  };  
  ///Set up the output stateSpace for this object
  ///
  virtual void addOptions(Options &opt,const string &prefix=""){
    GLens::addOptions(opt,prefix);
    addOption("remap_q","Use remapped mass-ratio coordinate.");
    addOption("q0","Prior max in q (with q>1) with remapped q0. Default=1e5/","1e5");
    addOption("GLB_gauss_q","Set to assume Gaussian (not flat) prior for log-q"); 
    addOption("GLB_rWide","Binary width/distance cuttoff for applying perturbed signle lens treatment (Einstein units). Default=5","5"); 
  };
  ///Set state parameters
  ///
  void setState(const state &st){
    //  logL separation (in log10 Einstein units)
    //  q mass ratio
    //  alignment angle phi0, (binary axis rel to trajectory direction) at closest approach point
    GLens::setState(st);
    //checkWorkingStateSpace();
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
  void remap_q(double q_ref_val=1e7){
    do_remap_q=true;
    q_ref=q_ref_val;
  };

};

#endif
