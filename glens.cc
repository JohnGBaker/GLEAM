//Gravitational lens equation for microlensing
//Written by John G Baker NASA-GSFC (2014-2016)

#include "glens.hh"
#include <gsl/gsl_poly.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include <algorithm>
#include <complex>
#include "omp.h"
#ifdef USE_KIND_16
#include <quadmath.h>
#else
#define __float128 long double
#endif

const bool fix_nr4_roots=true;
const bool inv_test_mode=false;
//const bool inv_test_mode=true;
extern bool debugint;
bool verbose=false;
bool test_result=false;
bool save_thetas_poly=true;
//bool save_thetas_poly=false;
//bool save_thetas_wide=true;  Strangely, this seems to provide no advantage.  Maybe a bug, but didn't see it. 
bool save_thetas_wide=false; 
//
// Interface to external Skowron&Gould fortran polynomial solver routine
//

#ifndef USE_KIND_16

extern "C" void cmplx_roots_5_(complex<double> roots[5], int &first_3_roots_order_changed, complex<double> poly[6], const int &polish_only);
void cmplx_roots_5(complex<double> roots[5], bool &first_3_roots_order_changed, complex<double> poly[6], const bool &polish_only){
  int first3=first_3_roots_order_changed;
  cmplx_roots_5_(roots, first3, poly, polish_only);
  first_3_roots_order_changed=first3;
};
extern "C" void cmplx_roots_gen_(complex<double> roots[], complex<double> poly[], const int &degree, const int &polish_roots_after, const int &use_roots_as_starting_points);
void cmplx_roots_gen(complex<double> roots[], complex<double> poly[], const int &degree, bool const &polish_roots_after, const bool &use_roots_as_starting_points){cmplx_roots_gen_(roots,poly,degree,polish_roots_after,use_roots_as_starting_points);};
typedef  long double ldouble;
const double LEADTOL=1e-5;
//const double LEADTOL=1e-4;
//const double epsTOL=1e-15;
const double epsTOL=1e-14;
//const double epsTOL=1e-11;

#else

extern "C" void cmplx_roots_5_(complex<__float128> roots[5], int &first_3_roots_order_changed, complex<__float128> poly[6], const int &polish_only);
void cmplx_roots_5(complex<double> roots[5], bool &first_3_roots_order_changed, complex<double> poly[6], const bool &polish_only){
  int first3=first_3_roots_order_changed;
  complex<__float128>longroots[5],longpoly[6];
  for(int i=0;i<6;i++)longpoly[i]=poly[i];
  if(!polish_only)for(int i=0;i<5;i++)longroots[i]=roots[i];
  //cout<<sizeof(longpoly)<<endl;
  cmplx_roots_5_(longroots, first3, longpoly, polish_only);
  for(int i=0;i<5;i++)roots[i]=longroots[i];  
  first_3_roots_order_changed=first3;
};
extern "C" void cmplx_roots_gen_(complex<__float128> roots[], complex<__float128> poly[], const int &degree, const int &polish_roots_after, const int &use_roots_as_starting_points);
void cmplx_roots_gen(complex<double> roots[], complex<double> poly[], const int &degree, const bool &polish_roots_after, const bool &use_roots_as_starting_points){
  complex<__float128>longroots[5],longpoly[6];
  for(int i=0;i<degree+1;i++)longpoly[i]=poly[i];
  /*char c[128];
  for(int i=0;i<degree+1;i++){
    quadmath_snprintf (c, sizeof c, "%-#*.35Qe", 46, longpoly[i]);  
    cout<<"poly["<<i<<"]="<<c<<endl;  
    }*/
  if(use_roots_as_starting_points)for(int i=0;i<degree;i++)longroots[i]=roots[i];
  //cout<<sizeof(longpoly)<<endl;
  int deg=degree;
  cmplx_roots_gen_(longroots,longpoly,deg,polish_roots_after,use_roots_as_starting_points);
  for(int i=0;i<degree;i++)roots[i]=longroots[i];  
};
typedef __float128 ldouble;
const double LEADTOL=3e-7;
//const double LEADTOL=1e-5;
const double epsTOL=1e-18;
//const double epsTOL=1e-14;
#endif



//
// ******************************************************************
// GLens routines ***************************************************
// ******************************************************************
//

//Compute curves through the 
void GLens::inv_map_curve(const vector<Point> &curve, vector<vector<Point> > &curve_images, vector<vector<double>> &curve_mags)
{
  // Given a curve of Points through the observer plane, compute the corresponding curve in the lens plane.
  //
  //Arguments:
  //
  //vector<Point> curve               -provides the list of observer-plane pointe
  //bool image_curves                 -the vectors comprising theta_series are each set to the same length and reordered
  //                                   so that matrix can be transposed to a matrix of curves this is achieved by finding
  //                                   nearest-neighbor alignments according to the ordering in curve, each theta will be
  //                                   of the same (maximum) length with parity=0 if there is no solution for that image
  //vector<vector<int>> & mag         -signed magnification for that image, or zero for no image
  //
  //

  //internal
  ///clear the outputs
  curve_images.clear();
  curve_mags.clear();

  //Main loop over curve
  int Ngrid=curve.size();
  for(int i=0; i<Ngrid;i++){
    Point beta=curve[i];
    vector<Point> thetas;
    //cout<<"curve: beta=("<<beta.x<<","<<beta.y<<")"<<endl;
    thetas.clear();
    thetas=invmap(beta);
    //record results;
    curve_images.push_back(thetas);
    vector<double> mags;
    for(Point th:thetas)mags.push_back(mag(th));
    curve_mags.push_back(mags);
  }
}

//Use GSL routine to integrate polygon trajectory
//just a sketch...
/*void GLens::integrate_invmap_curve (const vector<Point> &curve, vector<vector<Point> > &curve_images, vector<vector<double>> &curve_mags)
{
  //
  //control parameters:
  const double intTOL = GL_int_tol;  //control integration error tolerance
  const double test_result_tol = 1e-4;  //control integration error tolerance
  if(have_integrate)integrate=use_integrate;
  double prec=cout.precision();cout.precision(20);
  const double rWide_int_fac=100.0;

  ///clear the outputs
  time_series.clear();
  thetas_series.clear();
  index_series.clear();
  mag_series.clear();

  trajectory=&traj;//a convenience for passing to the integrator
  int NintSize=2*NimageMax;

  //Set up for GSL integration routines
  const gsl_odeiv2_step_type * stepType = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_step * step=nullptr;
  gsl_odeiv2_control * control=nullptr;
  gsl_odeiv2_evolve * evol=nullptr;
  gsl_odeiv2_system sys = {poly_root_integration_func_vec, NULL, (size_t)NintSize, this};
  double h=1e-5;
  step = gsl_odeiv2_step_alloc (stepType, NintSize);
  control = gsl_odeiv2_control_y_new (intTOL, 0.0);
  //control = gsl_odeiv2_control_standard_new (intTOL, 0.0,1.0,10.0);//Also apply limits on the derivative
  evol = gsl_odeiv2_evolve_alloc (NintSize);
  gsl_odeiv2_control_init(control,intTOL,0.0,1,0);//Do I need this?
  gsl_odeiv2_evolve_reset(evol);//or this
  
  
  //Main loop over observation times specified in the Trajectory object
  //initialization
  Point beta;
  vector<Point> thetas;
  double t_old;
  for(int i=0; i<Ngrid;i++){
    double tgrid=traj.get_obs_time(i);
    //entering main loop
    
    double t=t_old;
    //The next loop steps (probably more finely) toward the next grid time point.      
    while(t<tgrid){
      //integrate each image
      if(debugint)cout<<"\n t="<<t<<" mg="<<mg<<endl;
      Ntheta=thetas.size();
      double theta[NintSize];
      for(int image=0;image<thetas.size();image++){
	theta[2*image]=thetas[image].x;
	theta[2*image+1]=thetas[image].y;
      }
      for(int k=2*thetas.size();k<NintSize;k++)theta[k]=0;
      int status = gsl_odeiv2_evolve_apply (evol, control, step, &sys, &t, tgrid, &h, theta);
      //Need some step-size control checking for near caustics?
      if (status != GSL_SUCCESS) { 
	evolving=false;//switch to polynomial
	break;
      }
      //record results
      for(int image=0;image<thetas.size();image++){
	thetas[image]=Point(theta[2*image],theta[2*image+1]);	    
      }
      time_series.push_back(t);
      thetas_series.push_back(thetas);
    }

    gsl_odeiv2_evolve_reset(evol);
    
    
    if(time_series.size()<1)cout<<"Time series empty i="<<i<<endl;
    t_old=tgrid;
    index_series.push_back(time_series.size()-1);
  }//end of main observation times loop
  
  gsl_odeiv2_evolve_free (evol);
  gsl_odeiv2_control_free (control);
  gsl_odeiv2_step_free (step);
  
}
*/

///Assign a vector of subject points to the closest pairing of another vector of model points
///
///Result is returned as a vector of int of the same length as the subject vector, which must
///be no longer than the model vector.  When model is longer then a vector of the leftovers
//and a equal-length vector of the original index locations of these are return in leftovers
//and leftovers_map.
//In maxnorm, the function returns the largest squared-distance between any of the associated points.
vector<int> assign_points(const vector<Point> model, const vector<Point> subject, vector<Point> &leftovers,vector<int> &leftovers_map,double &maxnorm){
  //We assume the model contains at least as many elements as the subject so that all
  //subject element will be assigned to the nearest model element (or some overall nearest)
  //For now we do these in order.  It might be better to do a global fit...
  maxnorm=0;
  if(model.size()<subject.size()){
    cout<<"sizes: model="<<model.size()<<" subject="<<subject.size()<<endl;
    cout<<"(glens.cc)assign_points: model must include at least as many points as subject!"<<endl;
    exit(1);
  }
  vector<int> assignments(subject.size());
  vector<bool> used(model.size(),false);
  for(int i=0;i<subject.size();i++){
    double mindist2=1e100;
    int j0=-1;
    Point elem=subject[i];
    for(int j=0;j<model.size();j++){
      if(used[j])continue;
      Point cand=model[j];
      double dx=cand.x-elem.x;
      double dy=cand.y-elem.y;
      double dist2=dx*dx+dy*dy;
      if(j0<0||mindist2>dist2){
	mindist2=dist2;
	j0=j;
      }
      //cout<<"j="<<j<<" dist2="<<dist2<<" <= "<<mindist2<<" -> "<<j0<<endl;
    }
    if(mindist2>maxnorm)maxnorm=mindist2;
    assignments[i]=j0;
    used[j0]=true;
  }
  leftovers.clear();
  leftovers_map.clear();
  if(model.size()>subject.size()){
    for(int j=0;j<model.size();j++){
      if(not used[j]){
	leftovers.push_back(model[j]);
	leftovers_map.push_back(j);
      }
    }
  }
  //cout<<"returning:";for(int i=0;i<assignments.size();i++)cout<<" "<<assignments[i];cout<<endl;
  return assignments;  
}

///Some tools for the segment lists used in the image_area function below
typedef pair<int,int> joint;
bool jointOrder(const joint &j1, const joint &j2) {
  return j1.second < j2.second;
};

///Compute source image curves
///
///The source is described by a set of points defining a polygon.  The function maps that into a set of polygons defining
///images of the source. 
void GLens::compute_image_curves(const vector<Point> &polygon, const double maxlen, const double refine_limit, int & N, vector<vector<Point>> &closed_curves){
  ///polygon      : input  -closed polygon curve vertices
  ///radius       : input  -source radius sets scale for some of the precision cutoffs
  ///N            : output -Number of polygon points as refined
  //closed_curves : output -The image curves.

  vector<Point> curve=polygon;
  N=curve.size();

  double const twopi=2*M_PI;
  ///Controls
  const double expansion_limit=2.0;//we will refine if an image edge is more than this much times longer than the original polygon edge.
  const double maxnorm_limit=maxlen;
  const bool refine_sphere=false;
  const int nadd_max=3;

  //Get image points and mags
  vector<vector<Point> >  image_points;
  vector<vector<double> > image_point_mags;
  vector<double> vertex_mags;  //Used only for variance estimate at end.
  
  ///Step 1: Preliminary loop over the polygon vertices.
  ///This will be the most computationally intensive step unless we need to do a lot of refinement near caustics. For each source polygon vertex we compute a set of image points.
  ///
  for(int i=0; i<N;i++){
    Point beta=curve[i];
    vector<Point> thetas;
    thetas.clear();
    thetas=invmap(beta);
    vector<double> mags;
    for(Point th:thetas)mags.push_back(mag(th));
    double total_mg=0;
    for(double mg : mags){
      total_mg+=abs(mg);
    }
    //record results;
    image_points.push_back(thetas);
    image_point_mags.push_back(mags);
    vertex_mags.push_back(total_mg);
  }
    
  ///Step 2: We select a vertex to serve as a suitable starting place.
  ///        We need a starting vertex free from basic pathologies in its set of images.  In particular the number of images should be in the range between the minimum NimageMin and maximum NimageMax allowed values.  Images added beyond the minimum number should come in pairs (with opposite parities), so we require then n-NimageMin is even, and that the sum of parities has the expected value.  The set of image points from the selected vertex will provide the seeds from which we will construct (by tracing around the set of polygon vertices) a set of open image curves.  
  vector<vector<Point>>image_curves(NimageMax,vector<Point>(N+1));
  vector<vector<bool>>empty(NimageMax,vector<bool>(N+1,true));
  vector<vector<int>>mate(NimageMax,vector<int>(N+1,-1));
  vector<int>parities(NimageMax);
  int i0=-1;
  for(int i=0;i<N;i++){
    int ni=image_points[i].size();
    int netp=0;//parity sum
    for(int j=0;j<ni;j++){
      int p=copysign(1.0,image_point_mags[i][j]);
      if(image_point_mags[i][j]==0)p=-1;
      netp+=p;
    }
    vector<int>even_curves,odd_curves;
    if( ni >= NimageMin && ni <= NimageMax && (ni-NimageMin)%2==0 && -netp==NimageMin%2 ){
      i0=i;
      for(int j=0;j<NimageMax;j++){
	if(j<ni){ 
	  int p=copysign(1.0,image_point_mags[i0][j]);
	  if(image_point_mags[i0][j]==0)p=-1;//This will probably be a negative (near point lens) image.
	  image_curves[j][0]=image_points[i0][j];
	  parities[j]=p;
	  empty[j][0]=false;
	} else { //alternate parities for empty slots
	  parities[j]=(j%2)*2-1;
	  empty[j][0]=true;
	}
	if(parities[j]>0)even_curves.push_back(j);
	else odd_curves.push_back(j);
      }
      break;
    }
  }
  if(i0<0){//Fail if not found
    cout<<"finite_area_mag: Failing. No suitable starting point."<<endl;
    if(false){//Dump details for debugging
      cout<<"GLens::image_area_mag: Found no vertex with acceptable number of images."<<endl;
      cout<<"GLens::image_area_mag: Details:"<<endl;
      cout<<print_info()<<endl;
      cout<<"centers=";for(int k=1;k<NimageMin;k++)cout<<"  "<<getCenter(k).x<<" "<<getCenter(k).y;
      cout<<"\npoints[0]:"<<endl;
      int ithis=0;
      Point b=curve[ithis];
      cout<<"b="<<b.x<<" "<<b.y<<endl;
      for(int j=0;j<image_points[ithis].size();j++){
	Point bprime=map(image_points[ithis][j]);
	Point db=bprime-b;
	cout<<" "<<image_points[ithis][j].x<<" "<<image_points[j][ithis].y<<" "<<image_point_mags[ithis][j]
	    <<" bprime= "<<bprime.x<<" "<<bprime.y
	    <<" delta= "<<db.x<<" "<<db.y<<endl;
      }
    }
    N=0;
    return;
  }  
  
  ///Step 3: Construct open image curves.
  ///We step along the polygon vertices, checking the image sets, assigning them to open curves, refining the polygon edges as we go, as needed.
  ///
  ///As we step around the polygon, as long as the number of image points does not change, then the set of image points from the initial vertex should extend continuously into a set of image curves.  Each curve will be assembled consistently from either even or odd parity image points.
  ///
  ///As we consider each new vertex point, we first check whether the number of images has increased, decreased or remained the same. If it has remained the same, then we need only assign the unordered set of vertex images to the existing set of curves. We identify each new image point by even or odd parity, then, assign them one-by-one to the closest curve with the same parity.
  ///
  ///If the number of images increases, then, after assigning the closest points, there should be a matching number of points left over of each parity.  These will become the beginning of a additional new open image curves.
  ///
  ///If the number of images decreases, then there should be two of the image curves of opposite parity with ends very close together.  We terminate the two closest opposite-parity open image curves, and match the remaining points as in the other cases.  [Actually we haven't implemented it this way (though this way might be better.  Instead we proceed as before, associating points with nearest like-parity curves, and terminating those which are left over. We identify "mates" closest opposite-parity pairs among the terminating ends for reference later.]
  ///
  ///The actual process proceeds as follows: Begin with two sets, last_evens and last_odds, containing the images from the last sample point along the polygon.  Then find the images of the next polygon sample point and assign them to groups of even and odd parity.  Next perform basic sanity checks testing for the same conditions as require of the initial vertex images.  If an image point seems to be missing, then we try to add one, as done for the initial vertex, otherwise we try to fail gracefully by skipping ahead to the next polygon edge point.
  ///
  ///When attaching points to curves, we check that the distance of the jump from one image curve point to the next along that curve is not greater than expansion_limit (2 times the arclength of the source-edge arc segment represented by each polygon edge). If the image jump is larger, then we stop processing and instead refine the edge of the polygon by an integer factor as needed so that we linearly expect the condition to be met (unless this would exceed a refinment limit that we still need to explain.).
  ///
  ///The next step (3D) is to identify the mates of ends of image curve segments to beginnings of other image curve segments, where we understand the odd parity segments to run backward.  Need to elaborate...
  ///
  /// We can probably simplify this code (maybe a lot) by directly producing segments as we sort the point-images the first time.
  
  //initialization
  double norefine=false,refine_end=false;
  double maxnorm;
  vector<joint>segment_ends,segment_begins;
  vector<Point>leftover;
  vector<int>ileftover,imap;
  vector<Point>evens,odds,last_evens,last_odds,leftevens,leftodds;
  vector<int>last_evens_ind,last_odds_ind,evens_ind,odds_ind,ileftevens,ileftodds;
  //initialize evens/odds
  {
    int netp=0;
    for(int j=0;j<image_points[i0].size();j++){
      int p=copysign(1.0,image_point_mags[i0][j]);
      if(image_point_mags[i0][j]==0)p=-1;
      netp+=p;
      if(p>0){
	evens.push_back(image_points[i0][j]);
	evens_ind.push_back(j);
      } else {
	odds.push_back(image_points[i0][j]);
	odds_ind.push_back(j);
      }
    }
  }
  //Loop
  //Note: ithis labels the current index point in the image_points (or curves) 
  //      i labels the current index point in the image_curves
  int ilast=0,ithis=0;
  for(int i=1;i<=(image_points.size()) or refine_end;i++){
    //Step 3A: Initialization of loop interior
    //Note: By the end we will have gone through the loop once for each point that ends up in image_points
    //We may insert points into image_points
    //The starting point is i0 (Note we need to increment i0 if we insert points at/before i0;
    N=image_points.size();
    ilast=ithis;
    if( i==1 )ilast=(i+i0-1)%N;//FIXME Why doesn't this work???
    ithis=(i+i0)%N;
    //int last_ni=image_points[ilast].size();
    int last_ni=evens.size()+odds.size();//Unlike image_points[ilast].size() this is correct when the image checks fail for the last point.
    int ni=image_points[ithis].size();
    int netp=0;
    bool refine=false;
    
    //Step 3B: Assign points to even / odd lists
    last_evens=evens;last_evens_ind=evens_ind;evens.clear();evens_ind.clear();
    last_odds=odds;last_odds_ind=odds_ind;odds.clear();odds_ind.clear();
    netp=0;
    for(int j=0;j<ni;j++){
      int p=copysign(1.0,image_point_mags[ithis][j]);
      if(image_point_mags[ithis][j]==0)p=-1;
      netp+=p;
      if(p>0){
	evens.push_back(image_points[ithis][j]);
	evens_ind.push_back(j);
      } else {
	odds.push_back(image_points[ithis][j]);
	odds_ind.push_back(j);
      }
    }
    
    //special case: if we need to refine the end
    if(refine_end){
      refine_end=false;
    }

    //Step 3C: Basic sanity check on the next points image set
    bool pass=( ni >= NimageMin && ni <= NimageMax && (ni-NimageMin)%2==0 && -netp==NimageMin%2 );
    if(not pass){
      //Mitigation:
      //We will have to skip this point, but to minimize the damage, we first try to refine the previous edge (this will necessarily continue until max refinement is reached)
      if(not norefine){//can refine more
	refine=true;
      }else{
	//Can't refine more. Copy forward the last image point, and hope for the best...
	//for(int k=0;k<NimageMax;k++)image_curves[k][ithis]=image_curves[k][ilast];
	for(int k=0;k<NimageMax;k++)image_curves[k][i]=image_curves[k][i-1];
	evens=last_evens;evens_ind=last_evens_ind;
	odds=last_odds;odds_ind=last_odds_ind;
	norefine=false;
	continue;
      }
    } else { ///passed sanity check
      ///Step 3D: Assign new set of points to curves and check if refinement is needed 
      ///         We must handle connectivity given a variety of cases for the number of images.
      //         { same # of images, number increased, number decreased }
      for(int idummy=0;idummy<1;idummy++){//this is a dummy loop to enable conveniently break-ing out of the analysis if we need to refine
	if(ni<last_ni){
	  //Case: Number of images decreased.  First have to figure which to drop;
	  //idea: loop over all images at new point, and assign each to an old point then make those empty
	  //First we connect the even curves
	  if(last_evens.size()<evens.size())cout<<"Case 1a: last_evens<evens"<<endl;
	  imap=assign_points(last_evens,evens,leftevens,ileftevens,maxnorm);
	  if(maxnorm>maxnorm_limit and not norefine){refine=true;//cout<<"maxnorm 1a="<<maxnorm<<endl;
	    break;}
	  for(int k=0;k<evens.size();k++){
	    int ic=last_evens_ind[imap[k]];
	    image_curves[ic][i]=evens[k]; // add image to selected curve
	    empty[ic][i]=false;
	    evens_ind[k]=ic;
	  }
	  //Then the odd curves
	  imap=assign_points(last_odds,odds,leftodds,ileftodds,maxnorm);
	  if(maxnorm>maxnorm_limit and not norefine){refine=true;//cout<<"maxnorm 1b="<<maxnorm<<endl;
	    break;}
	  for(int k=0;k<odds.size();k++){
	    int ic=last_odds_ind[imap[k]];
	    image_curves[ic][i]=odds[k]; // add image to selected curve
	    empty[ic][i]=false;
	    odds_ind[k]=ic;
	  }	  
	  //Now mate the terminating ends
	  imap=assign_points(leftodds,leftevens,leftover,ileftover,maxnorm);
	  if(maxnorm>maxnorm_limit and not norefine){refine=true;//cout<<"maxnorm 1c="<<maxnorm<<endl;
	    break;}
	  for(int k=0;k<leftevens.size();k++){
	    int jeven=last_evens_ind[ileftevens[k]];
	    int jodd=last_odds_ind[ileftodds[imap[k]]];
	    mate[jeven][ilast]=jodd;
	    segment_ends.push_back(joint(jeven,(i-1)%N));
	    mate[jodd][ilast]=jeven;
	    segment_begins.push_back(joint(jodd,(i-1)%N));
	  }
	  
	} else {
	  // Case: Same number of images or new images have appeared;
	  imap=assign_points(evens,last_evens,leftevens,ileftevens,maxnorm);
	  if(maxnorm>maxnorm_limit and not norefine){refine=true;//cout<<"maxnorm 2a="<<maxnorm<<endl;
	    break;}
	  for(int k=0;k<last_evens.size();k++){
	    int ic=last_evens_ind[k];
	    image_curves[ic][i]=evens[imap[k]]; // add image to selected curve
	    empty[ic][i]=false;
	    evens_ind[imap[k]]=ic;
	  }
	  imap=assign_points(odds,last_odds,leftodds,ileftodds,maxnorm);
	  if(maxnorm>maxnorm_limit and not norefine){refine=true;//cout<<"maxnorm 2b="<<maxnorm<<endl;
	    break;}
	  for(int k=0;k<last_odds.size();k++){
	    int ic=last_odds_ind[k];
	    image_curves[ic][i]=odds[imap[k]]; // add image to selected curve
	    empty[ic][i]=false;
	    odds_ind[imap[k]]=ic;
	  }
	  //now mate the terminating ends  (Note: at this point it should be true that leftevens and leftodds are equal-length
	  if(leftodds.size()!=leftevens.size()){
	    cout<<"This shouldn't happen. New odds not equal to new evens!"<<endl;
	    exit(1);
	  }
	  if(leftodds.size()>0){
	    //Case : New images have appeared
	    imap=assign_points(leftodds,leftevens,leftover,ileftover,maxnorm);
	    if(maxnorm>maxnorm_limit and not norefine){refine=true;
	      break;}
	    int jeven=0,jodd=0;
	    for(int k=0;k<leftevens.size();k++){
	      for(int j=jeven;j<NimageMax;j++){
		jeven=j;
		if(empty[j][i] and parities[j]>0)break;
	      }
	      image_curves[jeven][i]=leftevens[k];
	      empty[jeven][i]=false;
	      for(int j=jodd;j<NimageMax;j++){
		jodd=j;
		if(empty[j][i] and parities[j]<0)break;
	      }
	      image_curves[jodd][i]=leftodds[k];
	      empty[jodd][i]=false;
	      mate[jeven][i]=jodd;
	      segment_begins.push_back(joint(jeven,i));
	      mate[jodd][i]=jeven;
	      segment_ends.push_back(joint(jodd,i));
	      evens_ind[ileftevens[k]]=jeven;
	      odds_ind[ileftodds[k]]=jodd;
	    }
	  }
	}
      }
    }

    //Step 3E: Refine if needed, as determined in above process
    if(refine){
      //We refine by integer an integer factor
      //After refinement want: maxnorm -> maxnorm/factor^2 < maxnorm_limit
      //So we need: factor^2 > maxnorm/maxnorm_limit
      int nadd=sqrt(maxnorm/maxnorm_limit);
      if(nadd>nadd_max)nadd=nadd_max;//enforce a maximum degree of refinement in one step
      double factor=1+nadd;
      Point p0=curve[ilast];
      Point dp=curve[ithis]-curve[ilast];
      double dp2=dp.x*dp.x+dp.y*dp.y;
      while(nadd>0 and dp2<refine_limit*factor*factor){nadd--;factor--;}//enforce an overall limit in how far we refine
      if(nadd==0)norefine=true;//We will go through this step again with a flag set indicating not to refine more
      double phi0=0,dphi=0;
      //We assemble everthing we need into vectors, then insert these into the originals
      vector<vector<Point> > new_image_points;
      vector<vector<double> >new_image_point_mags;
      vector<Point> new_curve_points;
      int iinsert=ithis;
      if(iinsert==0)iinsert=N;
      for(int k=0;k<nadd;k++){
	Point pnew;
	new_curve_points.push_back(pnew);
	vector<Point> thetas=invmap(pnew);
	new_image_points.push_back(thetas);
	vector<double> mags;
	for(Point th:thetas)mags.push_back(mag(th));
	new_image_point_mags.push_back(mags);
      }
      curve.insert(curve.begin()+iinsert,new_curve_points.begin(),new_curve_points.end());
      image_points.insert(image_points.begin()+iinsert,new_image_points.begin(),new_image_points.end());
      image_point_mags.insert(image_point_mags.begin()+iinsert,new_image_point_mags.begin(),new_image_point_mags.end());
      //also need to add space in the image curves for the new points.
      for(int j=0;j<NimageMax;j++){
	image_curves[j].resize(N+nadd+1);
	empty[j].resize(N+nadd+1,true);
	mate[j].resize(N+nadd+1,-1);
      }
      if(i0>=iinsert)i0+=nadd;
      
      //Having executed the refinement, now step back to retry (new) next segment
      i--;
      evens=last_evens;evens_ind=last_evens_ind;
      odds=last_odds;odds_ind=last_odds_ind;
      ithis=ilast;
    }
    else norefine=false;  //reset after each step
    
    //Step 3F TBD???
    //(This is an experimental effort to address the "should be adding a point" message issue
    //but it needs work, causes immediate crash now.)
    //This is an early start on the preparation for step 4, which requires
    //tieing the starts and ends together.  
    if(false and i==(image_points.size())){ //We only need to do this if we have reached the loop end.
      vector<Point>end,start;
      vector<int>iend,istart;
      for(int j=0;j<NimageMax;j++){
	if(not empty[j][0]){
	  start.push_back(image_curves[j][0]);
	  istart.push_back(j);
	}
	if(not empty[j][N]){
	  end.push_back(image_curves[j][N]);
	  iend.push_back(j);
	}
      } 
      if(start.size()<end.size())cout<<"start<end"<<endl;
      imap=assign_points(start,end,leftover,ileftover,maxnorm);
      if(maxnorm>maxnorm_limit){
	cout<<"We *are* adding a point!"<<endl;
	//we need to add it in the next loop cycle when the start point is seen as the "next" one
	if(not norefine)refine_end=true;
      }
    }
  }//end of step 3 loop constructing image curves
  N=image_points.size();

  //Step 4: Define segment connections to sew the ends together.
  //Now we define segment connections to the ends together
  //We have looped to N above, one extra for closure, but we match these to the 0 slot.
  //These will be the same set of points as at index 0, but may not be assigned to the same curves;
  vector<Point>end,start;
  vector<int>iend,istart;
  for(int j=0;j<NimageMax;j++){
    if(not empty[j][0]){
      start.push_back(image_curves[j][0]);
      istart.push_back(j);
    }
    if(not empty[j][N]){
      end.push_back(image_curves[j][N]);
      iend.push_back(j);
    }
  }
  //cout<<"end:";for(int i=0;i<end.size();i++)cout<<"  "<<end[i].x<<" "<<end[i].y;cout<<endl;
  //cout<<"start:";for(int i=0;i<start.size();i++)cout<<"  "<<start[i].x<<" "<<start[i].y;cout<<endl;
  if(start.size()<end.size())cout<<"start<end"<<endl;
  imap=assign_points(start,end,leftover,ileftover,maxnorm);
  if(maxnorm>maxnorm_limit){
    cout<<"We should be adding a point!"<<endl;
  }
  for(int j=0;j<end.size();j++){
    int jend=iend[j];
    int jstart=istart[imap[j]];
    if(mate[jend][N]<0){//no mate already assigned at end
      if(parities[jend]>0){
	segment_ends.push_back(joint(jend,N-1));
	segment_begins.push_back(joint(jstart,0));
      } else {
	segment_ends.push_back(joint(jstart,0));
	segment_begins.push_back(joint(jend,N-1));
      }
    } else {
      //Special case that we have identified the N-label point as a start.
      //This should really be at 0...
      //but we must it may be on  a different curve's 0, so we have to fix it now.
      //cout<<"mate"<<endl;
      for(int k=0;k<segment_ends.size();k++){
	if(segment_ends[k]==joint(jend,N)){
	  segment_ends[k]=joint(jstart,0);
	}
	if(segment_begins[k]==joint(jend,N)){
	  segment_begins[k]=joint(jstart,0);
	}
      }
    }
  }
  
  // Step 5:
  //Next we connect the segments (abstractly) into some number of closed curves as vectors of segment-joint indices
  vector<vector<int>>closed_curve_segments;//
  vector<bool>used(segment_ends.size(),false);
  for(int i=0;i<segment_begins.size();i++){//loop over all the 
    if(used[i])continue;
    //start a new curve, beginning with this segment
    closed_curve_segments.push_back(vector<int>(0));
    vector<int>this_curve;
    this_curve.push_back(i);
    used[i]=true;
    int seg=i;
    while(true){
      int jbegin=segment_begins[seg].first;
      int ibegin=segment_begins[seg].second;
      //find where this segment ends:
      seg=-1;
      int mindelta=N+1;
      for(int k=0;k<segment_ends.size();k++){
	if(segment_ends[k].first!=jbegin)continue;//j must match
	int iend=segment_ends[k].second;
	int delta=(iend-ibegin)*parities[jbegin];
	if(delta<0)continue;//must end after starting (in parity dirirection)
	if(delta<mindelta){
	  mindelta=delta;
	  seg=k;
	}
      }
      if(seg<0){
	if(this_curve.size()<2)
	  cout<<"Failed to find second point in curve"<<endl;
	else
	  cout<<"Failed to find next segment. Something is wrong!:"<<endl;	
      }
      if(seg<0 or used[seg]){
	closed_curve_segments.back()=this_curve;
	break;//this segment already listed, closes off curve (other ways to test, but this seems fast).
      } else {
	this_curve.push_back(seg);
	used[seg]=true;
      }
    }
  }

  // Step 6: Now we actually assemble the closed curves
  closed_curves.clear();
  for(auto scurve:closed_curve_segments){
    vector<Point>curve;
    for(int i=0;i<scurve.size();i++){
      int iseg=scurve[i];
      int inext=scurve[(i+1)%scurve.size()];
      int j=segment_begins[iseg].first;
      int istart=segment_begins[iseg].second;
      int iend=segment_ends[inext].second;
      if(parities[j]>0){
	//(Seems this should work for C++11 standard even for neg parity, but with gcc5 seems not.) 
	curve.insert(curve.end(),image_curves[j].begin()+istart,image_curves[j].begin()+iend+parities[j]);
      } else {
	//neg parity
	//first insert in the wrong order
	curve.insert(curve.end(),image_curves[j].begin()+iend,image_curves[j].begin()+istart+1);
	//then reverse order in place
	reverse(curve.end()-istart+iend-1,curve.end());
      }
    }
    closed_curves.push_back(curve);
  }  
}

///Compute magnification of a finite N-sided polygon by finding the area of images based on the method Gould-Gaucherel
///
///The method requires the following steps: First make a list of points defining the polygon, an apply inv_map_curve to
///find the associated image points.  Next we compute the infinitessimal (signed) magnification of point.  Next we must
///organize these points into curves.  There are a couple general principles to this.  First, each image has a parity
///sign, which should be the same as the sign of each of its component point infinitessimal magnifications.  Second,
///though there should always be an even number of image points beyond NimageMin, when the polygon spans a caustic, some
///images may not exist for all polygon points.  If the latter is not true then we try to infill a missing point or
///remove a spurious image point. When there is a jump in the number of images then the added or subtracted images must
///be paired.  The association of curves into points is realized by a minimum-area principle with the constraint that
///parities must align.
///
///The Gould-Goucherel method recognizes the image curves as either closed curves, or paired segments.  The points are
///regarded as vectors from some origin. For each edge within a closed curve or interior section of segment p2-p1 , the area
///contribution is computed by \f$dA=p1*p2\f$ (or \f$dA=p2*p1\f$ for negative parity).  Though the individual segments contributions will
///depend on the choice of origin, the overall area will not.  Note that when the origin is outside the image, the extra area
///will be compensated by negative area somewhere else.  The edges of the segments are must be sewn together by finding pairs
///or opposing-parity curves with an edge at the same vertex.  These curves are, in effect, stiched together by adding area
///contributions as would arise if the curve was connected \f$dA=p1_{pos} \times p2_{neg}\f$.  The magnification is computed by dividing the
///summed image area by the similarly computed polygon area.
///
///The function returns an estimate for the 'variance' over the surface.
///Presently, this is simply the varance of the point-magnifications around the polygon
///wrt the area magnification.
///
void GLens::image_area_mag(Point &p, double radius, int & N, double &magnification, double &var, ostream *out,vector<vector<Point> > *outcurves){
  ///p      : input  -observer position, relative source center
  ///       : output -returns centroid shift
  ///radius : input  -source radius
  ///N      : input  -Number of polygon points
  ///       : output -Number of sampled points
  ///magnification   : output magnification result
  ///var    : output estimated variance result
  ///out    : input  : optional output stream for dumping curve results
  ///outcurves: output: optional return of the computed closed curve set.
  //save these for debugging reference 
  Point p0=p;
  int N0=N;

  double const twopi=2*M_PI;

  ///Controls
  //const double expansion_limit=1.05;//we will refine if an image edge is more than this much times longer than the original polygon edge.
  const double expansion_limit=2.0;//we will refine if an image edge is more than this much times longer than the original polygon edge.
  //const double expansion_limit=10.0;//we will refine if an image edge is more than this much times longer than the original polygon edge.
  const double maxnorm_limit=finite_source_tol+pow(expansion_limit*2*M_PI*radius/N,2.0);
  //const double maxnorm_limit=finite_source_tol;
  const int nadd_max=3;
  const bool refine_sphere=false;
  const double refine_prec_limit=1e-12;
  //const double refine_prec_limit=1e-15;
  const double refine_limit=pow(twopi*radius/N/finite_source_refine_limit,2.0)+(p.x*p.x+p.y*p.y)*refine_prec_limit*refine_prec_limit;
  bool debug_area_mag=false;
  bool super_debug=false;
  bool fix_vertex_images=true;
  //N*=20;

  ///special debugging

  //double dq=183261.3878-dynamic_cast<GLensBinary*>(this)->get_q();
  //double dL=3.48323-dynamic_cast<GLensBinary*>(this)->get_L();
  //double px0=1.36179711562218;
  //double px0=-1.454926;
  //#pragma omp critical	
  //{
  //cout<<"p.x="<<p.x<<" "<<print_info()<<endl;
  //cout<<"px0="<<px0<<" dq="<<dq<<" dL="<<dL<<endl;
  //}
  //if(false and abs(dq)<.001 and abs(dL) and abs(p.x-px0)<.001){
  //cout<<"Special debugging"<<endl;
  //debug_area_mag=true;}
  
  //Construct polygon
  double dphi=twopi/N;
  vector<Point> curve(N);
  for(int i=0;i<N;i++)curve[i]=Point(p.x+radius*cos(dphi*i),p.y+radius*sin(dphi*i));
  
  //Get image points and mags
  vector<vector<Point> >  image_points;
  vector<vector<double> > image_point_mags;
  vector<double> vertex_mags;  //Used only for variance estimate at end.
  //inv_map_curve(curve,image_points,image_point_mags);
  
  ///Step 1: Preliminary loop over the polygon vertices.
  ///This will be the most computationally intensive step unless we need to do a lot of refinement near caustics. For each source polygon vertex we compute a set of image points.
  ///
  ///Then (unless fix_vertex_image==false) we check for errors in the image set and try to fix them. In particular, if there are fewer than NimageMin (eg 3 for binarly lens) image points, then we try to recover one.  We are most likely to have missed a negative parity image point which is extremely close to one of the lens centers.  If so, the sum of parities will be zero.  If it is, then we determine which lens center is farthest from any of the image points, and add one there. (How often does this occur?)
  for(int i=0; i<N;i++){
    Point beta=curve[i];
    vector<Point> thetas;
    thetas.clear();
    thetas=invmap(beta);
    vector<double> mags;
    for(Point th:thetas)mags.push_back(mag(th));
    double total_mg=0;
    for(double mg : mags){
      total_mg+=abs(mg);
    }
    if(fix_vertex_images){ //check for some kinds of errors and try to fix:
      int ni=thetas.size();
      if( NimageMin==ni+1 ){//Seem to be missing an image
	//Do some prep
	int netp=0;
	vector<Point>odds;
	for(int j=0;j<thetas.size();j++){
	  int p=copysign(1.0,mags[j]);
	  if(mags[j]==0)p=-1;
	  netp+=p;
	  if(p<0)odds.push_back(thetas[j]);
	}
	if(netp==0){//Seems we can fix it!
	  //Try to detect which lens center is farthest from any of the odd images:
	  int kmiss=-1;
	  double rknormfar=0;
	  for(int k=1;k<NimageMin;k++){
	    Point c=getCenter(k);
	    //find closest odd image to this center
	    double rknormmin=INFINITY;
	    for(auto b:odds){
	      Point rk=b+c*(-1.0);
	      double rknorm=rk.x*rk.x+rk.y*rk.y;
	      if(rknorm<rknormmin)rknormmin=rknorm;
	    }
	    if(rknormfar<rknormmin){
	      rknormfar=rknormmin;
	      kmiss=k;
	    }
	  }
	  Point cmiss=getCenter(kmiss);
	  cout<<"VERTEX lens_center "<<kmiss<<" at ("<<cmiss.x<<","<<cmiss.y<<") seems to be missing its image; adding an image at this center!"<<endl;	
	  thetas.push_back(cmiss);
	  mags.push_back(-1e-100);
	}
      }
    }
    //record results;
    image_points.push_back(thetas);
    image_point_mags.push_back(mags);
    vertex_mags.push_back(total_mg);
  }
  /*
  cout<<"first pass image points:"<<endl;
  for(int j=0;j<N;j++){
    cout<<j<<" ";
    for(int k=0;k<image_points[j].size();k++){
      cout<<"  "<<image_points[j][k].x<<" "<<image_points[j][k].y;
    }
    cout<<endl;
  }
  */
    
  ///Step 2: We select a vertex to serve as a suitable starting place.
  ///        We need a starting vertex free from basic pathologies in its set of images.  In particular the number of images should be in the range between the minimum NimageMin and maximum NimageMax allowed values.  Images added beyond the minimum number should come in pairs (with opposite parities), so we require then n-NimageMin is even, and that the sum of parities has the expected value.  The set of image points from the selected vertex will provide the seeds from which we will construct (by tracing around the set of polygon vertices) a set of open image curves.  
  vector<vector<Point>>image_curves(NimageMax,vector<Point>(N+1));
  vector<vector<bool>>empty(NimageMax,vector<bool>(N+1,true));
  vector<vector<int>>mate(NimageMax,vector<int>(N+1,-1));
  vector<int>parities(NimageMax);
  int i0=-1;
  for(int i=0;i<N;i++){
    int ni=image_points[i].size();
    int netp=0;//parity sum
    for(int j=0;j<ni;j++){
      int p=copysign(1.0,image_point_mags[i][j]);
      if(image_point_mags[i][j]==0)p=-1;
      netp+=p;
    }
    //cout<<i<<" < "<<N<<" ni="<<ni<<" NimageMin="<<NimageMin<<" NimageMax="<<NimageMax<<" netp="<<netp<<endl;
    //cout<<" pass=( "<<(ni >= NimageMin)<<" && "<<(ni <= NimageMax)<<" && "<<((ni-NimageMin)%2==0)<<" && "<<(-netp==NimageMin%2)<<" )"<<endl;
    vector<int>even_curves,odd_curves;
    if( ni >= NimageMin && ni <= NimageMax && (ni-NimageMin)%2==0 && -netp==NimageMin%2 ){
      i0=i;
      for(int j=0;j<NimageMax;j++){
	if(j<ni){ 
	  int p=copysign(1.0,image_point_mags[i0][j]);
	  if(image_point_mags[i0][j]==0)p=-1;//This will probably be a negative (near point lens) image.
	  image_curves[j][0]=image_points[i0][j];
	  parities[j]=p;
	  empty[j][0]=false;
	} else { //alternate parities for empty slots
	  parities[j]=(j%2)*2-1;
	  //cout<<"set parities["<<j<<"]="<<parities[j]<<" j%2="<<(j%2)<<endl;
	  empty[j][0]=true;
	}
	if(parities[j]>0)even_curves.push_back(j);
	else odd_curves.push_back(j);
      }
      break;
    }
    //cout<<"not pass"<<endl;
  }
  if(i0<0){//Fail if not found
    magnification=INFINITY;
    var=0;
    cout<<"finite_area_mag: Failing. No suitable starting point."<<endl;
    if(false){//Dump details for debugging
      cout<<"GLens::image_area_mag: Found no vertex with acceptable number of images."<<endl;
      cout<<"GLens::image_area_mag: Details:"<<endl;
      cout<<print_info()<<endl;
      cout<<"centers=";for(int k=1;k<NimageMin;k++)cout<<"  "<<getCenter(k).x<<" "<<getCenter(k).y;
      cout<<"\npoints[0]:"<<endl;
      int ithis=0;
      Point b=curve[ithis];
      cout<<"b="<<b.x<<" "<<b.y<<endl;
      for(int j=0;j<image_points[ithis].size();j++){
	Point bprime=map(image_points[ithis][j]);
	Point db=bprime-b;
	cout<<" "<<image_points[ithis][j].x<<" "<<image_points[j][ithis].y<<" "<<image_point_mags[ithis][j]
	    <<" bprime= "<<bprime.x<<" "<<bprime.y
	    <<" delta= "<<db.x<<" "<<db.y<<endl;
      }
    }
    return;
  }  
  //static int i0save=-1;if(i0!=i0save)cout<<" i0="<<i0<<endl;i0save=i0;
  if(debug_area_mag)cout<<"Enter: i0="<<i0<<endl;
  //cout<<" parities: "<<parities[0]<<" "<<parities[1]<<" "<<parities[2]<<" "<<parities[3]<<" "<<parities[4]<<endl;

  
  ///Step 3: Construct open image curves.
  ///We step along the polygon vertices, checking the image sets, assigning them to open curves, refining the polygon edges as we go, as needed.
  ///
  ///As we step around the polygon, as long as the number of image points does not change, then the set of image points from the initial vertex should extend continuously into a set of image curves.  Each curve will be assembled consistently from either even or odd parity image points.
  ///
  ///As we consider each new vertex point, we first check whether the number of images has increased, decreased or remained the same. If it has remained the same, then we need only assign the unordered set of vertex images to the existing set of curves. We identify each new image point by even or odd parity, then, assign them one-by-one to the closest curve with the same parity.
  ///
  ///If the number of images increases, then, after assigning the closest points, there should be a matching number of points left over of each parity.  These will become the beginning of a additional new open image curves.
  ///
  ///If the number of images decreases, then there should be two of the image curves of opposite parity with ends very close together.  We terminate the two closest opposite-parity open image curves, and match the remaining points as in the other cases.  [Actually we haven't implemented it this way (though this way might be better.  Instead we proceed as before, associating points with nearest like-parity curves, and terminating those which are left over. We identify "mates" closest opposite-parity pairs among the terminating ends for reference later.]
  ///
  ///The actual process proceeds as follows: Begin with two sets, last_evens and last_odds, containing the images from the last sample point along the polygon.  Then find the images of the next polygon sample point and assign them to groups of even and odd parity.  Next perform basic sanity checks testing for the same conditions as require of the initial vertex images.  If an image point seems to be missing, then we try to add one, as done for the initial vertex, otherwise we try to fail gracefully by skipping ahead to the next polygon edge point.
  ///
  ///When attaching points to curves, we check that the distance of the jump from one image curve point to the next along that curve is not greater than expansion_limit (2 times the arclength of the source-edge arc segment represented by each polygon edge). If the image jump is larger, then we stop processing and instead refine the edge of the polygon by an integer factor as needed so that we linearly expect the condition to be met (unless this would exceed a refinment limit that we still need to explain.).
  ///
  ///The next step (3D) is to identify the mates of ends of image curve segments to beginnings of other image curve segments, where we understand the odd parity segments to run backward.  Need to elaborate...
  
  //initialization
  double norefine=false,refine_end=false;
  double maxnorm;
  vector<joint>segment_ends,segment_begins;
  vector<Point>leftover;
  vector<int>ileftover,imap;
  vector<Point>evens,odds,last_evens,last_odds,leftevens,leftodds;
  vector<int>last_evens_ind,last_odds_ind,evens_ind,odds_ind,ileftevens,ileftodds;
  //initialize evens/odds
  {
    int netp=0;
    for(int j=0;j<image_points[i0].size();j++){
      int p=copysign(1.0,image_point_mags[i0][j]);
      if(image_point_mags[i0][j]==0)p=-1;
      netp+=p;
      if(p>0){
	evens.push_back(image_points[i0][j]);
	evens_ind.push_back(j);
      } else {
	odds.push_back(image_points[i0][j]);
	odds_ind.push_back(j);
      }
    }
  }
  //Loop
  //Note: ithis labels the current index point in the image_points (or curves) 
  //      i labels the current index point in the image_curves
  int ilast=0,ithis=0;
  for(int i=1;i<=(image_points.size()) or refine_end;i++){
    //Step 3A: Initialization of loop interior
    //Note: By the end we will have gone through the loop once for each point that ends up in image_points
    //We may insert points into image_points
    //The starting point is i0 (Note we need to increment i0 if we insert points at/before i0;
    N=image_points.size();
    ilast=ithis;
    if( i==1 )ilast=(i+i0-1)%N;//FIXME Why doesn't this work???
    ithis=(i+i0)%N;
    //cout<<"i ilast ithis:"<<i<<" "<<ilast<<" "<<ithis<<endl;
    //int last_ni=image_points[ilast].size();
    int last_ni=evens.size()+odds.size();//Unlike image_points[ilast].size() this is correct when the image checks fail for the last point.
    int ni=image_points[ithis].size();
    int netp=0;
    bool refine=false;
    
    //Step 3B: Assign points to even / odd lists
    last_evens=evens;last_evens_ind=evens_ind;evens.clear();evens_ind.clear();
    last_odds=odds;last_odds_ind=odds_ind;odds.clear();odds_ind.clear();
    netp=0;
    for(int j=0;j<ni;j++){
      int p=copysign(1.0,image_point_mags[ithis][j]);
      if(image_point_mags[ithis][j]==0)p=-1;
      netp+=p;
      if(p>0){
	evens.push_back(image_points[ithis][j]);
	evens_ind.push_back(j);
      } else {
	odds.push_back(image_points[ithis][j]);
	odds_ind.push_back(j);
      }
    }
    
    //special case: if we need to refine the end
    if(refine_end){
      refine_end=false;
    }

    //Step 3C: Basic sanity check on the next points image set
    if(debug_area_mag)cout<<"i="<<i<<" ilast="<<ilast<<" ithis="<<ithis<<" i0="<<i0<<" norefine="<<norefine<<endl;
    bool pass=( ni >= NimageMin && ni <= NimageMax && (ni-NimageMin)%2==0 && -netp==NimageMin%2 );
    //cout<<" ni="<<ni<<"  last_ni="<<last_ni<<endl;
    //cout<<"i="<<i<<" ni="<<ni<<" NimageMin="<<NimageMin<<" NimageMax="<<NimageMax<<" netp="<<netp<<endl;
    //cout<<"  pass=( "<<(ni >= NimageMin)<<" && "<<(ni <= NimageMax)<<" && "<<((ni-NimageMin)%2==0)<<" && "<<(-netp==NimageMin%2)<<" )"<<endl;
    if(not pass){
      //cout<<"not pass."<<endl;
      if(NimageMin-ni==1 and netp==0){//Detect lost an image point near one of the point lenses
	//Try to detect which lens center is farthest from any of the odd images:
	int kmiss=-1;
	double rknormfar=0;
	for(int k=1;k<NimageMin;k++){
	  Point c=getCenter(k);
	  //find closest odd image to this center
	  double rknormmin=INFINITY;
	  for(auto b:odds){
	    Point rk=b+c*(-1.0);
	    double rknorm=rk.x*rk.x+rk.y*rk.y;
	    if(rknorm<rknormmin)rknormmin=rknorm;
	  }
	  if(rknormfar<rknormmin){
	    rknormfar=rknormmin;
	    kmiss=k;
	  }
	}
	Point cmiss=getCenter(kmiss);
	//cout<<"lens_center "<<kmiss<<" at ("<<cmiss.x<<","<<cmiss.y<<") seems to be missing its image; adding an image at this center!"<<endl;	
	//cout<<"N="<<N<<endl;
	image_points[ithis].push_back(cmiss);
	image_point_mags[ithis].push_back(-1e-100);
	odds.push_back(cmiss);
	odds_ind.push_back(ni);
	ni++;
	netp=-1;
      } else { //Did not pass check and did not manage an easy fix.
	if(debug_area_mag){
	  cout<<"GLens::image_area_mag: Point did not pass need to add handling of unexpected patterns of images."<<endl;
	  cout<<"N="<<N<<endl;
	  cout<<"ni="<<ni<<" NimageMin="<<NimageMin<<" NimageMax="<<NimageMax<<" netp="<<netp<<endl;
	  cout<<" pass=( "<<(ni >= NimageMin)<<" && "<<(ni <= NimageMax)<<" && "<<((ni-NimageMin)%2==0)<<" && "<<(-netp==NimageMin%2)<<" )"<<endl;
	  cout<<print_info()<<endl;
	  cout<<"centers=";for(int k=1;k<NimageMin;k++)cout<<"  "<<getCenter(k).x<<" "<<getCenter(k).y;
	  cout<<"\npoints:"<<endl;
	  Point b=curve[ithis];
	  cout<<"b="<<b.x<<" "<<b.y<<endl;
	  for(int j=0;j<image_points[ithis].size();j++){
	    Point db=map(image_points[ithis][j])-b;
	    cout<<" "<<image_points[ithis][j].x<<" "<<image_points[ithis][j].y<<" "<<image_point_mags[ithis][j]<<" delta= "<<db.x<<" "<<db.y<<endl;
	  }
	  cout<<"last points:"<<endl;
	  cout<<" ilast="<<ilast<<" < "<<image_points.size()<<endl;
	  cout<<"b="<<b.x<<" "<<b.y<<endl;
	  for(int j=0;j<image_points[ilast].size();j++){
	    Point db=map(image_points[ilast][j])-b;
	    cout<<" "<<image_points[ilast][j].x<<" "<<image_points[ilast][j].y<<" "<<image_point_mags[ilast][j]<<" delta= "<<db.x<<" "<<db.y<<endl;
	  }
	  cout<<"last_evens"<<endl;for(int j=0;j<last_evens.size();j++)cout<<last_evens[j].x<<" "<<last_evens[j].y<<" "<<mag(last_evens[j])<<endl;
	  cout<<"last_odds"<<endl;for(int j=0;j<last_odds.size();j++)cout<<last_odds[j].x<<" "<<last_odds[j].y<<" "<<mag(last_odds[j])<<endl;
	  cout<<"  ...trying to fail gracefully by neglecting this point..."<<endl;
	}
	//debug_area_mag=true;
	//Mitigation:
	//We will have to skip this point, but to minimize the damage, we first try to refine the previous edge (this will necessarily continue until max refinement is reached)
	if(not norefine){//can refine more
	  refine=true;
	}else{
	  //Can't refine more. Copy forward the last image point, and hope for the best...
	  //for(int k=0;k<NimageMax;k++)image_curves[k][ithis]=image_curves[k][ilast];
	  for(int k=0;k<NimageMax;k++)image_curves[k][i]=image_curves[k][i-1];
	  evens=last_evens;evens_ind=last_evens_ind;
	  odds=last_odds;odds_ind=last_odds_ind;
	  norefine=false;
	  continue;
	}
      }
    }

    if(not refine){
      ///Step 3D: Assign new set of points to curves and check if refinement is needed 
      ///         We must handle connectivity given a variety of cases for the number of images.
      //         { same # of images, number increased, number decreased }
      if(debug_area_mag){
	cout<<"#odds/even sizes: "<<odds.size()<<"/"<<evens.size()<<"   images="<<image_points[ithis].size();
	cout<<"  "<<empty[0][ilast]<<empty[1][ilast]<<empty[2][ilast]<<empty[3][ilast]<<empty[4][ilast]<<endl;
	cout<<" last_evens_ind("<<last_evens.size()<<"): ";for(int j=0;j<last_evens_ind.size();j++)cout<<" "<<last_evens_ind[j];cout<<endl;
	cout<<" last_odds_ind("<<last_odds.size()<<"):  ";for(int j=0;j<last_odds_ind.size();j++)cout<<" "<<last_odds_ind[j];cout<<endl;
      }
      for(int idummy=0;idummy<1;idummy++){//this is a dummy loop to enable conveniently break-ing out of the analysis if we need to refine
	if(ni<last_ni){
	  //Case 1: number of images decreased.  First have to figure which to drop;
	  //idea: loop over all images at new point, and assign each to an old point then make those empty
	  //First we connect the even curves
	  if(debug_area_mag)cout<<"case1a"<<endl;
	  if(last_evens.size()<evens.size())cout<<"Case 1a: last_evens<evens"<<endl;
	  imap=assign_points(last_evens,evens,leftevens,ileftevens,maxnorm);
	  if(maxnorm>maxnorm_limit and not norefine){refine=true;//cout<<"maxnorm 1a="<<maxnorm<<endl;
	    break;}
	  for(int k=0;k<evens.size();k++){
	    int ic=last_evens_ind[imap[k]];
	    //cout<<"evens["<<k<<"] -> image_curves["<<ic<<"]"<<endl;
	    image_curves[ic][i]=evens[k]; // add image to selected curve
	    empty[ic][i]=false;
	    evens_ind[k]=ic;
	  }
	  //Then the odd curves
	  if(debug_area_mag)cout<<"case1b"<<endl;
	  if(last_odds.size()<odds.size())cout<<"Case 1b: last_odds<odds"<<endl;
	  imap=assign_points(last_odds,odds,leftodds,ileftodds,maxnorm);
	  if(maxnorm>maxnorm_limit and not norefine){refine=true;//cout<<"maxnorm 1b="<<maxnorm<<endl;
	    break;}
	  for(int k=0;k<odds.size();k++){
	    int ic=last_odds_ind[imap[k]];
	    //cout<<"odds["<<k<<"] -> image_curves["<<ic<<"]"<<endl;
	    image_curves[ic][i]=odds[k]; // add image to selected curve
	    empty[ic][i]=false;
	    odds_ind[k]=ic;
	  }	  
	  //Now mate the terminating ends
	  if(debug_area_mag)cout<<"case1c"<<endl;
	  if(leftodds.size()<leftevens.size())cout<<"Case 1c: model<subject"<<endl;
	  imap=assign_points(leftodds,leftevens,leftover,ileftover,maxnorm);
	  if(maxnorm>maxnorm_limit and not norefine){refine=true;//cout<<"maxnorm 1c="<<maxnorm<<endl;
	    break;}
	  for(int k=0;k<leftevens.size();k++){
	    int jeven=last_evens_ind[ileftevens[k]];
	    int jodd=last_odds_ind[ileftodds[imap[k]]];
	    //cout<<"leftevens["<<k<<"] -> image_curves["<<jeven<<"]"<<endl;
	    //cout<<"leftodds["<<k<<"] -> image_curves["<<jodd<<"]"<<endl;
	    mate[jeven][ilast]=jodd;
	    segment_ends.push_back(joint(jeven,(i-1)%N));
	    mate[jodd][ilast]=jeven;
	    segment_begins.push_back(joint(jodd,(i-1)%N));
	    ////this is wrong...
	    //evens_ind[ileftevens[k]]=jeven;
	    //odds_ind[ileftodds[k]]=jodd;
	  }
	  
	} else {
	  // Cases 2/3: same number of images or new images have appeared;
	  //cout<<"Assign:"<<endl;
	  if(debug_area_mag)cout<<"case2a"<<endl;
	  if(evens.size()<last_evens.size()){
	    cout<<"Case 2a: model<subject"<<endl;
	    cout<<"p.x="<<p.x<<endl;
	    cout<<" pass=( "<<(ni >= NimageMin)<<" && "<<(ni <= NimageMax)<<" && "<<((ni-NimageMin)%2==0)<<" && "<<(-netp==NimageMin%2)<<" )"<<endl;
	    cout<<print_info()<<endl;
	    cout<<"centers=";for(int k=1;k<NimageMin;k++)cout<<"  "<<getCenter(k).x<<" "<<getCenter(k).y;
	    cout<<"\npoints:"<<endl;
	    Point b=curve[ithis];
	    cout<<"b="<<b.x<<" "<<b.y<<endl;
	    for(int j=0;j<image_points[ithis].size();j++){
	      Point db=map(image_points[ithis][j])-b;
	      cout<<" "<<image_points[ithis][j].x<<" "<<image_points[ithis][j].y<<" "<<image_point_mags[ithis][j]<<" delta= "<<db.x<<" "<<db.y<<endl;
	    }
	    cout<<"evens"<<endl;for(int j=0;j<evens.size();j++)cout<<evens[j].x<<" "<<evens[j].y<<" "<<mag(evens[j])<<endl;
	    cout<<"odds"<<endl;for(int j=0;j<odds.size();j++)cout<<odds[j].x<<" "<<odds[j].y<<" "<<mag(odds[j])<<endl;
	    cout<<"last points:"<<endl;
	    cout<<"i="<<i<<" ilast="<<ilast<<" < "<<image_points.size()<<endl;
	    //b=curve[ilast];
	    cout<<"b="<<b.x<<" "<<b.y<<endl;
	    for(int j=0;j<image_points[ilast].size();j++){
	      Point db=map(image_points[ilast][j])-b;
	      cout<<" "<<image_points[ilast][j].x<<" "<<image_points[ilast][j].y<<" "<<image_point_mags[ilast][j]<<" delta= "<<db.x<<" "<<db.y<<endl;
	    }
	    cout<<"last_evens"<<endl;for(int j=0;j<last_evens.size();j++)cout<<last_evens[j].x<<" "<<last_evens[j].y<<" "<<mag(last_evens[j])<<endl;
	    cout<<"last_odds"<<endl;for(int j=0;j<last_odds.size();j++)cout<<last_odds[j].x<<" "<<last_odds[j].y<<" "<<mag(last_odds[j])<<endl;
	    cout<<print_info(20)<<endl;
	    cout.precision(20);
	    cout<<"Entered image_area_mag with: p0=("<<p0.x<<","<<p0.y<<")\n radius="<<radius<<" N0="<<N0<<endl;
	  }
	  imap=assign_points(evens,last_evens,leftevens,ileftevens,maxnorm);
	  if(maxnorm>maxnorm_limit and not norefine){refine=true;//cout<<"maxnorm 2a="<<maxnorm<<endl;
	    break;}
	  for(int k=0;k<last_evens.size();k++){
	    int ic=last_evens_ind[k];
	    //cout<<"evens["<<imap[k]<<"] -> image_curves["<<ic<<"]"<<endl;
	    image_curves[ic][i]=evens[imap[k]]; // add image to selected curve
	    empty[ic][i]=false;
	    evens_ind[imap[k]]=ic;
	  }
	  if(debug_area_mag)cout<<"case2b"<<endl;
	  if(odds.size()<last_odds.size())cout<<"Case 2b: model<subject"<<endl;
	  imap=assign_points(odds,last_odds,leftodds,ileftodds,maxnorm);
	  if(maxnorm>maxnorm_limit and not norefine){refine=true;//cout<<"maxnorm 2b="<<maxnorm<<endl;
	    break;}
	  for(int k=0;k<last_odds.size();k++){
	    int ic=last_odds_ind[k];
	    //cout<<"odds["<<imap[k]<<"] -> image_curves["<<ic<<"]"<<endl;
	    image_curves[ic][i]=odds[imap[k]]; // add image to selected curve
	    empty[ic][i]=false;
	    odds_ind[imap[k]]=ic;
	  }
	  //now mate the terminating ends  (Note: at this point it should be true that leftevens and leftodds are equal-length
	  if(leftodds.size()!=leftevens.size()){
	    cout<<"This shouldn't happen. New odds not equal to new evens!"<<endl;
	    exit(1);
	  }
	  if(leftodds.size()>0){
	    //Case 3: New images have appeared
	    //cout<<"New!"<<endl;
	    //cout<<"#odds/even sizes: "<<odds.size()<<"/"<<evens.size()<<"   images="<<image_points[ithis].size();
	    //cout<<"  "<<empty[0][ilast]<<empty[1][ilast]<<empty[2][ilast]<<empty[3][ilast]<<empty[4][ilast]<<endl;
	    //cout<<" last_odds_ind("<<last_odds.size()<<"):  ";for(int j=0;j<last_odds_ind.size();j++)cout<<" "<<last_odds_ind[j];cout<<endl;
	    //cout<<" last_evens_ind("<<last_evens.size()<<"):  ";for(int j=0;j<last_evens_ind.size();j++)cout<<" "<<last_evens_ind[j];cout<<endl;
	    if(debug_area_mag)cout<<"case3"<<endl;
	    if(leftodds.size()<leftevens.size())cout<<"Case 3: model<subject"<<endl;
	    imap=assign_points(leftodds,leftevens,leftover,ileftover,maxnorm);
	    if(maxnorm>maxnorm_limit and not norefine){refine=true;//cout<<"maxnorm 3="<<maxnorm<<endl;
	      break;}
	    int jeven=0,jodd=0;
	    for(int k=0;k<leftevens.size();k++){
	      for(int j=jeven;j<NimageMax;j++){
		jeven=j;
		if(empty[j][i] and parities[j]>0)break;
	      }
	      //cout<<"leftevens["<<k<<"] -> image_curves["<<jeven<<"]"<<endl;
	      image_curves[jeven][i]=leftevens[k];
	      empty[jeven][i]=false;
	      for(int j=jodd;j<NimageMax;j++){
		jodd=j;
		if(empty[j][i] and parities[j]<0)break;
	      }
	      //cout<<"leftodds["<<k<<"] -> image_curves["<<jodd<<"]"<<endl;
	      image_curves[jodd][i]=leftodds[k];
	      empty[jodd][i]=false;
	      mate[jeven][i]=jodd;
	      segment_begins.push_back(joint(jeven,i));
	      mate[jodd][i]=jeven;
	      segment_ends.push_back(joint(jodd,i));
	      evens_ind[ileftevens[k]]=jeven;
	      odds_ind[ileftodds[k]]=jodd;
	    }
	  }
	}
      }
    }

    //Step 3E: Refine if called for
    if(refine){
      //We refine by integer an integer factor
      //After refinement want: maxnorm -> maxnorm/factor^2 < maxnorm_limit
      //So we need: factor^2 > maxnorm/maxnorm_limit
      int nadd=sqrt(maxnorm/maxnorm_limit);
      if(nadd>nadd_max)nadd=nadd_max;//enforce a maximum degree of refinement in one step
      Point p0=curve[ilast];
      Point dp=curve[ithis]-curve[ilast];
      double dp2=dp.x*dp.x+dp.y*dp.y;
      double factor=1+nadd;
      while(nadd>0 and dp2<refine_limit*factor*factor){nadd--;factor--;}//enforce an overall limit in how far we refine
      if(nadd==0)norefine=true;//We will go through this step again with a flag set indicating not to refine more
      if(debug_area_mag){
	cout<<(norefine?"NOT ":"")<<"Refining! N="<<N<<" -> "<<N+nadd<<" maxnorm="<<maxnorm<<" > "<<maxnorm_limit<<endl;
       
	cout<<"f,dp2,refine_limit,f^2:"<<factor<<","<<dp2<<","<<refine_limit<<","<<factor*factor<<endl;
	cout<<"p=("<<p.x<<","<<p.y<<")  radius="<<radius<<endl;
	
	cout<<"p(last)=("<<curve[ilast].x<<","<<curve[ilast].y<<")"<<endl;
	cout<<"p(this)=("<<curve[ithis].x<<","<<curve[ithis].y<<")"<<endl;
	cout<<"Before:"<<endl;
	cout<<"evens"<<endl;for(int j=0;j<evens.size();j++)cout<<evens[j].x<<" "<<evens[j].y<<" "<<mag(evens[j])<<endl;
	cout<<"odds"<<endl;for(int j=0;j<odds.size();j++)cout<<odds[j].x<<" "<<odds[j].y<<" "<<mag(odds[j])<<endl;
	cout<<"last_evens"<<endl;for(int j=0;j<last_evens.size();j++)cout<<last_evens[j].x<<" "<<last_evens[j].y<<" "<<mag(last_evens[j])<<endl;
	cout<<"last_odds"<<endl;for(int j=0;j<last_odds.size();j++)cout<<last_odds[j].x<<" "<<last_odds[j].y<<" "<<mag(last_odds[j])<<endl;
      }
      double phi0=0,dphi=0;
      if(refine_sphere){
	dp=p0-p;
	phi0=atan2(dp.y,dp.x);
	dp=curve[ithis]-p;
	double phi1=atan2(dp.y,dp.x);
	if(phi1<phi0)phi1+=2*M_PI;
	dphi-phi1-phi0;
      }
      //We assemble everthing we need into vectors, then insert these into the originals
      vector<vector<Point> > new_image_points;
      vector<vector<double> >new_image_point_mags;
      vector<Point> new_curve_points;
      int iinsert=ithis;
      if(iinsert==0)iinsert=N;
      for(int k=0;k<nadd;k++){
	Point pnew;
	if(refine_sphere){
	  double phinew= phi0 + dphi * ((k+1.0)/factor);
	  pnew=p+Point(cos(phinew),sin(phinew))*radius;
	} else pnew = p0 + dp * ((k+1.0)/factor);
	new_curve_points.push_back(pnew);
	vector<Point> thetas=invmap(pnew);
	new_image_points.push_back(thetas);
	vector<double> mags;
	for(Point th:thetas)mags.push_back(mag(th));
	new_image_point_mags.push_back(mags);
	//cout<<"inserting point ("<<pnew.x<<","<<pnew.y<<") at position "<<iinsert+k<<endl;
	//cout<<" thetas = ";for(int j=0;j<thetas.size();j++)cout<<thetas[j].x<<" "<<thetas[j].y<<"  ";cout<<endl;
      }
      /*
      cout<<"before inserting:"<<endl;
      for(int j=0;j<N;j++){
	cout<<j<<" ";
	for(int k=0;k<image_points[j].size();k++){
	  cout<<"  "<<image_points[j][k].x<<" "<<image_points[j][k].y;
	}
	cout<<endl;
	}
      */
      curve.insert(curve.begin()+iinsert,new_curve_points.begin(),new_curve_points.end());
      image_points.insert(image_points.begin()+iinsert,new_image_points.begin(),new_image_points.end());
      image_point_mags.insert(image_point_mags.begin()+iinsert,new_image_point_mags.begin(),new_image_point_mags.end());
      //also need to add space in the image curves for the new points.
      for(int j=0;j<NimageMax;j++){
	image_curves[j].resize(N+nadd+1);
	empty[j].resize(N+nadd+1,true);
	mate[j].resize(N+nadd+1,-1);
      }
      if(i0>=iinsert)i0+=nadd;
      
      if(debug_area_mag){
	cout<<"inserted "<<nadd<<" points. i0="<<i0<<" len(image_points)="<<image_points.size()<<endl;
	cout<<"After:"<<endl;
	cout<<"evens"<<endl;for(int j=0;j<evens.size();j++)cout<<evens[j].x<<" "<<evens[j].y<<" "<<mag(evens[j])<<endl;
	cout<<"odds"<<endl;for(int j=0;j<odds.size();j++)cout<<odds[j].x<<" "<<odds[j].y<<" "<<mag(odds[j])<<endl;
	cout<<"last_evens"<<endl;for(int j=0;j<last_evens.size();j++)cout<<last_evens[j].x<<" "<<last_evens[j].y<<" "<<mag(last_evens[j])<<endl;
	cout<<"last_odds"<<endl;for(int j=0;j<last_odds.size();j++)cout<<last_odds[j].x<<" "<<last_odds[j].y<<" "<<mag(last_odds[j])<<endl;
	/*
	cout<<"after inserting:"<<endl;
	for(int j=0;j<N+nadd;j++){
	  cout<<j<<" ";
	  for(int k=0;k<image_points[j].size();k++){
	    cout<<"  "<<image_points[j][k].x<<" "<<image_points[j][k].y;
	  }
	  cout<<endl;
	}
	*/
      }
      
      //Now we try to forget this step ever happened
      i--;
      evens=last_evens;evens_ind=last_evens_ind;
      odds=last_odds;odds_ind=last_odds_ind;
      ithis=ilast;
    }
    else norefine=false;  //reset after each step
    if(debug_area_mag){
      cout<<"done with i="<<i<<" segment_ends.size="<<segment_ends.size()<<endl<<endl;
      cout<<"end of loop interior."<<endl; 
      /*
	for(int j=0;j<N;j++){
	cout<<j<<" ";
	for(int k=0;k<image_curves.size();k++){
	  cout<<"  "<<image_curves[k][j].x<<" "<<image_curves[k][j].y;
	  if(empty[k][j])cout<<"[empty]";
	}
	cout<<endl;
      }
      */
    }

    //Step 3F 
    //(This is an experimental effort to address the "should be adding a point" message issue
    //but it needs work, causes immediate crash now.)
    //This is an early start on the preparation for step 4, which requires
    //tieing the starts and ends together.  
    if(false and i==(image_points.size())){ //We only need to do this if we have reached the loop end.
      vector<Point>end,start;
      vector<int>iend,istart;
      for(int j=0;j<NimageMax;j++){
	if(not empty[j][0]){
	  start.push_back(image_curves[j][0]);
	  istart.push_back(j);
	}
	if(not empty[j][N]){
	  end.push_back(image_curves[j][N]);
	  iend.push_back(j);
	}
      }
      if(start.size()<end.size())cout<<"start<end"<<endl;
      imap=assign_points(start,end,leftover,ileftover,maxnorm);
      if(maxnorm>maxnorm_limit){
	cout<<"We *are* adding a point!"<<endl;
	//we need to add it in the next loop cycle when the start point is seen as the "next" one
	debug_area_mag=true;
	if(not norefine)refine_end=true;
      }
    }
  }//end of step 3 loop constructing image curves
  N=image_points.size();

  if(super_debug){
    cout<<" N --> "<<N<<endl;
    cout<<"final image points:"<<endl;
    for(int j=0;j<N;j++){
      cout<<j<<" ";
      for(int k=0;k<image_points[j].size();k++){
	cout<<"  "<<image_points[j][k].x<<" "<<image_points[j][k].y;
      }
      cout<<endl;
    }
    cout<<"initial image curves:"<<endl;
    for(int j=0;j<N;j++){
      cout<<j<<" ";
      for(int k=0;k<image_curves.size();k++){
	cout<<"  "<<image_curves[k][j].x<<" "<<image_curves[k][j].y;
	if(empty[k][j])cout<<"[empty]";
      }
      cout<<endl;
    }
    cout<<" parities:";
    for(int i=0;i<parities.size();i++)cout<<"  "<<parities[i]<<endl;
    cout<<"Segments:"<<endl;
    for(int k=0;k<segment_ends.size();k++){
      cout<<"  end:"<<k<<"=("<<segment_ends[k].first<<","<<segment_ends[k].second<<")"<<endl;
      cout<<"  begin:"<<k<<"=("<<segment_begins[k].first<<","<<segment_begins[k].second<<")"<<endl;
    }
  }

  //Step 4: Define segment connections to sew the ends together.
  //Now we define segment connections to the ends together
  //We have looped to N above, one extra for closure, but we match these to the 0 slot.
  //These will be the same set of points as at index 0, but may not be assigned to the same curves;
  vector<Point>end,start;
  vector<int>iend,istart;
  for(int j=0;j<NimageMax;j++){
    if(not empty[j][0]){
      start.push_back(image_curves[j][0]);
      istart.push_back(j);
    }
    if(not empty[j][N]){
      end.push_back(image_curves[j][N]);
      iend.push_back(j);
    }
  }
  //cout<<"end:";for(int i=0;i<end.size();i++)cout<<"  "<<end[i].x<<" "<<end[i].y;cout<<endl;
  //cout<<"start:";for(int i=0;i<start.size();i++)cout<<"  "<<start[i].x<<" "<<start[i].y;cout<<endl;
  if(start.size()<end.size())cout<<"start<end"<<endl;
  imap=assign_points(start,end,leftover,ileftover,maxnorm);
  if(maxnorm>maxnorm_limit){
    cout<<"We should be adding a point!"<<endl;
    debug_area_mag=true;
  }
  for(int j=0;j<end.size();j++){
    int jend=iend[j];
    int jstart=istart[imap[j]];
    //cout<<"j="<<j<<" jend="<<jend<<" jstart="<<jstart<<endl;
    if(mate[jend][N]<0){//no mate already assigned at end
      //cout<<"no mate"<<endl;
      if(parities[jend]>0){
	segment_ends.push_back(joint(jend,N-1));
	segment_begins.push_back(joint(jstart,0));
      } else {
	segment_ends.push_back(joint(jstart,0));
	segment_begins.push_back(joint(jend,N-1));
      }
    } else {
      //Special case that we have identified the N-label point as a start.
      //This should really be at 0...
      //but we must it may be on  a different curve's 0, so we have to fix it now.
      //cout<<"mate"<<endl;
      for(int k=0;k<segment_ends.size();k++){
	//cout<<"k="<<k<<": "<<(segment_ends[k]==joint(jend,N))<<"  "<<(segment_begins[k]==joint(jend,N))<<endl;
	if(segment_ends[k]==joint(jend,N)){
	  segment_ends[k]=joint(jstart,0);
	}
	if(segment_begins[k]==joint(jend,N)){
	  segment_begins[k]=joint(jstart,0);
	}
      }
    }
  }
  //cout<<"Fixed ends"<<endl; 
  //for(int k=0;k<segment_ends.size();k++){
  //  cout<<"  end:"<<k<<"=("<<segment_ends[k].first<<(parities[segment_ends[k].first]>0?"+":"-")<<","<<segment_ends[k].second<<")"<<endl;
  // cout<<"  begin:"<<k<<"=("<<segment_begins[k].first<<(parities[segment_begins[k].first]>0?"+":"-")<<","<<segment_begins[k].second<<")"<<endl;
  // }
  
  // Step 5:
  //Next we connect the segments (abstractly) into some number of closed curves as vectors of segment-joint indices
  vector<vector<int>>closed_curve_segments;//
  vector<bool>used(segment_ends.size(),false);
  for(int i=0;i<segment_begins.size();i++){//loop over all the 
    if(used[i])continue;
    //cout<<"starting curve with seg "<<i<<endl;
    //start a new curve, beginning with this segment
    closed_curve_segments.push_back(vector<int>(0));
    vector<int>this_curve;
    this_curve.push_back(i);
    used[i]=true;
    int seg=i;
    while(true){
      int jbegin=segment_begins[seg].first;
      int ibegin=segment_begins[seg].second;
      //cout<<"seg begins:("<<jbegin<<","<<ibegin<<")"<<endl;
      //find where this segment ends:
      seg=-1;
      int mindelta=N+1;
      for(int k=0;k<segment_ends.size();k++){
	//cout<<"ends/begins:"<<segment_ends.size()<<"/"<<segment_begins.size()<<endl;
	//cout<<"trying seg:"<<k<<"=("<<segment_ends[k].first<<","<<segment_ends[k].second<<")"<<endl;
	if(segment_ends[k].first!=jbegin)continue;//j must match
	int iend=segment_ends[k].second;
	int delta=(iend-ibegin)*parities[jbegin];
	//cout<<"  ...delta="<<delta<<endl;
	if(delta<0)continue;//must end after starting (in parity dirirection)
	if(delta<mindelta){
	  //cout<<"  ...<mindelta="<<mindelta<<endl;
	  mindelta=delta;
	  seg=k;
	}
      }
      if(seg<0){
	if(this_curve.size()<2)
	  cout<<"Failed to find second point in curve"<<endl;
	else
#pragma omp critical	
	  {
	    cout<<"Failed to find next segment. Something is wrong!:"<<endl;
	    cout<<"p.x="<<p.x<<" "<<print_info()<<endl;
	    for(int k=0;k<segment_ends.size();k++){
	      cout<<k<<"b=("<<segment_ends[k].first<<","<<segment_ends[k].second<<")"<<endl;
	      cout<<k<<"e=("<<segment_begins[k].first<<","<<segment_begins[k].second<<")"<<endl;
	    }
	    for(int j=0;j<N;j++){
	      cout<<j<<" ";
	      for(int k=0;k<image_curves.size();k++){
		cout<<"  "<<image_curves[k][j].x<<" "<<image_curves[k][j].y;
		if(empty[k][j])cout<<"[empty]";
	      }
	      cout<<endl;
	    }
	  }
      }
      //cout<<" found next seg "<<seg<<endl;
      if(seg<0 or used[seg]){
	//cout<<"completed curve: ";
	/*for(int seg: this_curve){
	  cout<<"("<<segment_ends[seg].first<<","<<segment_ends[seg].second<<")="
	    <<"("<<segment_begins[seg].first<<","<<segment_begins[seg].second<<") ... ";
	}
	cout<<endl;*/
	closed_curve_segments.back()=this_curve;
	break;//this segment already listed, closes off curve (other ways to test, but this seems fast).
      } else {
	//cout<<" appending next seg "<<seg<<endl;
	this_curve.push_back(seg);
	used[seg]=true;
      }
      //cout<<"end of inner while loop: seg="<<seg<<endl;
    }
  }
  if(super_debug){
    cout<<"Computed curves"<<endl; 
    cout<<"# Found "<<closed_curve_segments.size()<<" closed curves:"<<endl;
    for(auto curve: closed_curve_segments){
      cout<<"#";
      for(auto seg: curve){
	cout<<"("<<segment_ends[seg].first<<","<<segment_ends[seg].second<<(parities[segment_ends[seg].first]>0?"+":"-")<<")="
	    <<"("<<segment_begins[seg].first<<","<<segment_begins[seg].second<<(parities[segment_begins[seg].first]>0?"+":"-")<<") ... ";
      }
      cout<<endl;
    }
  }

  // Step 6: Now we actually assemble the closed curves
  vector<vector<Point>>closed_curves;
  for(auto scurve:closed_curve_segments){
    vector<Point>curve;
    for(int i=0;i<scurve.size();i++){
      int iseg=scurve[i];
      int inext=scurve[(i+1)%scurve.size()];
      int j=segment_begins[iseg].first;
      int istart=segment_begins[iseg].second;
      int iend=segment_ends[inext].second;
      //cout<<"i,j,iseg,inext,istart,iend: "<<i<<" "<<j<<" "<<iseg<<" "<<inext<<" "<<istart<<" "<<iend<<endl;
      //cout<<"("<<segment_ends[iseg].first<<","<<segment_ends[iseg].second<<(parities[segment_ends[iseg].first]>0?"+":"-")<<")="
      //  <<"("<<segment_begins[iseg].first<<","<<segment_begins[iseg].second<<(parities[segment_begins[iseg].first]>0?"+":"-")<<") ... ";
      if(parities[j]>0){
	//cout<<"parities[j]="<<parities[j]<<" image_curves[j].size()="<<image_curves[j].size()<<endl;
	//(Seems this should work for C++11 standard even for neg parity, but with gcc5 seems not.) 
	curve.insert(curve.end(),image_curves[j].begin()+istart,image_curves[j].begin()+iend+parities[j]);
	//cout<<"inserted"<<endl;
      } else {
	//neg parity
	//first insert in the wrong order
	//cout<<"inserting"<<endl;
	//cout<<" image_curves[j].size()="<<image_curves[j].size()<<endl;
	curve.insert(curve.end(),image_curves[j].begin()+iend,image_curves[j].begin()+istart+1);
	//cout<<"reversing"<<endl;
	//then reverse order in place
	reverse(curve.end()-istart+iend-1,curve.end());
      }
    }
    //cout<<"curve "<<closed_curves.size()<<" is "<<endl;
    //for(int k=0;k<curve.size();k++)cout<<k<<" "<<curve[k].x<<" "<<curve[k].y<<endl;
    closed_curves.push_back(curve);
  }

  //copy curves to output vector if provided
  if(outcurves)*outcurves=closed_curves;
  
  //cout<<"Assembled curves"<<endl;
  if(out){
    *out<<endl;
    int n=0;
    for(auto & curve : closed_curves)if(curve.size()>n)n=curve.size();
    for(int i=0;i<n+1;i++){			
      for(auto & curve : closed_curves){
	int ii=i%curve.size();//just keep wrapping if curve is shorter than longest
	*out<<curve[ii].x<<" "<<curve[ii].y<<" ";
      }
      *out<<endl;
    }
    *out<<endl;
  }

  
  // Step 7: Compute area and overall image centroid
  double area=0;
  Point cent(0,0);
  for(int i=0;i<closed_curves.size();i++){
    vector<Point>&curve=closed_curves[i];
    Point pcentroid(0,0);//pointwise centroid
    for(auto pt :curve){
      pcentroid = pcentroid + pt;
    }
    pcentroid = pcentroid*(1.0/curve.size());
    //cout<<"pcentroid["<<i<<"]=("<<pcentroid.x<<","<<pcentroid.y<<")"<<endl;
    Point dAreaCentroid;
    double carea=getPolygonAreaCoM(curve,dAreaCentroid,pcentroid);
    //cout<<"curve area="<<area<<endl;
    //cout<<i<<" darea_centroid = "<<dAreaCentroid.x<<" "<<dAreaCentroid.y<<endl;
    area+=carea;
    cent=cent+dAreaCentroid;
    //this next part is a quality control diagnostic, to verify that we don't run into surprises where the signed area is not as expected...
    if(false){
      bool inside=false;
      for(int j=0;j<closed_curves.size();j++){
	if(j!=i and pointInPolygon(curve[0],closed_curves[j])!=0){
	  if(carea>0){
	    cout<<"curve "<<i<<" is inside curve "<<j<<" area["<<i<<"]="<<carea<<endl;
	    cout<<"GLens::image_area_mag: Area seems to have wrong sign!"<<endl;
	  } else inside=true;
	}
      }
      if(carea<0 and not inside){
	cout<<"curve "<<i<<" is NOT inside another curve but area["<<i<<"]="<<carea<<endl;
	cout<<"GLens::image_area_mag: Area seems to have wrong sign!"<<endl;
	cout<<"winding numbers for each curve around first point in polygon["<<i<<"]"<<endl;
	for(int j=0;j<closed_curves.size();j++)
	  cout<<pointInPolygon(curve[0],closed_curves[j])<<endl;
      }
      cout<<"count_area="<<carea<<endl;
    }  //End of diagnositic

  }
  Point c0;
  double area0=getPolygonAreaCoM(curve,p,c0);
  magnification=area/area0;
  p=cent*(1.0/area)+c0*(-1.0/area0);//we are recycling this to use as the overall image centroid offset now.
  //cout<<"total area="<<area<<" > "<<area0<<endl;
  //cout<<"centroid shift=("<<p.x<<","<<p.y<<")"<<endl;

  // Step 8: Compute variance, trimming excesses
  var=0;
  for( auto mg : vertex_mags){
    double dmg=mg-magnification;
    if(dmg>magnification)dmg=magnification;
    dmg*=source_var;
    var+=dmg*dmg;
  }
  var/=vertex_mags.size();
  //cout<<"var="<<var<<" n="<<vertex_mags.size()<<endl;
  //cout<<"sqrt(var)="<<sqrt(var)<<endl;
  
}


///Brute force integration finite size magnification around 1-D circle
///
///This is used below in computing the 2D integral over the image plane.
int GLens::brute_force_circle_mag(const Point &p, const double radius, const double ctol, double &magnification){
  //const double tol=finite_source_tol;
  const double twopi=2*M_PI;
  const double tol=ctol;
  const double tolcutmin=1e-12;
  //const double tolcutmin=1e-8;
  const double dphimin=twopi*tolcutmin/1e20;
  //const double tolcutmin=1e-12;
  //const double dphimin=twopi*tol;
  //int nphi=4*(2+1/sqrt(tol));
  int nphimin=20;
  int nphi=(nphimin+1/sqrt(tol));
  //int nphi=3;
  int nrefine=2;
  const double maxmag=1e15;
  
  //cout<<"circle r="<<radius<<endl;

  //set up the initial grid in rho
  vector<double> mags(nphi,0);
  vector<double> phis(nphi);
  for( int i=0;i<nphi;i++)phis[i]=i*2.0*M_PI/nphi;

  double oldnphi=0;
  //main loop of refinement passes over the set of radii
  magnification=1;
  int met_tol_count=0;
  while(oldnphi<nphi){
    oldnphi=nphi;
    //compute necessary circle magnifications
    for( int i=0;i<nphi;i++){
      if(mags[i]==0){
	double x=cos(phis[i])*radius;
	double y=sin(phis[i])*radius;
	Point beta=p+Point(x,y);
	vector<Point> thetas=invmap(beta);
	double intens=1.0;//can make this a function of r,phi for general intensity profile
	double mg=mag(thetas);
	if(not (mg<maxmag))mg=maxmag;
	mags[i]=mg*intens;
      }
    }

    if(1){
      //Now compute the mean mag
      //cout<<"PHIS ";
      //for(auto ph : phis)cout<<ph<<" ";
      //cout<<endl;
      double msum=0,psum=0;
      for(int i=0;i<nphi;i++){
	int il=((i-1) % nphi + nphi) % nphi;
	double dphi=phis[i]-phis[il];
	if(dphi<0)dphi+=twopi;
	if(dphi>twopi)dphi-=twopi;
	psum+=dphi;
	//cout<<il<<" "<<i<<" "<<dphi<<" "<<mags[i]<<" "<<mags[il]<<endl;
	msum+=(mags[i]+mags[il])*dphi/2.0;
      }
      //cout<<"psum="<<psum<<endl;
      msum/=twopi;
      //cout<<"    mag "<<msum<<" dmag="<<msum-magnification<<endl;    
      //if(msum-magnification>1){//diagnostic
      //cout<<"JUMP: phis= ";
      //for(auto ph : phis)cout<<ph<<" ";
      //cout<<endl;
      //}
      if(fabs(msum-magnification)<tol){
	if(met_tol_count>1)break;
	met_tol_count++;
      } else met_tol_count=0;
      magnification=msum;
    }
    
    vector<double> newmags;
    vector<double> newphis;
    //double margin=sqrt(nphi);
    //double margin=nphi;
    double margin=30;
    double tolcut=tol*twopi/margin;
    if(tolcut<tolcutmin)tolcut=tolcutmin;
    //initialize
    int i0=nphi-2;if(i0<0)i0+=nphi;
    double mlast=mags[i0],phlast=phis[i0];
    double mthis=mags[nphi-1],phthis=phis[nphi-1];
    bool norefine=false;
    bool norefinelast=false;
    
    for(int i=0;i<nphi;i++){
      if(i+1>=nphi){if(norefinelast)norefine=true;}
      //For the error estimate, Dm, the integral will be exact if the magnitude trends linearly
      //so the estimate Dm is equal to difference from the the linear interpolant
      //The contribution to the total is also ultimately proportional to dphi [how we improve near a discont.]
      double mnext=mags[i],phnext=phis[i];
      if(phthis>phnext)phthis-=twopi;
      if(phlast>phthis)phlast-=twopi;
      double dphi=phnext-phlast;
      double Dm=(dphi*mthis-(-phlast*mnext+phnext*mlast + phthis*(mnext-mlast)))/2.0;
      //cout<<i<<" "<<phlast<<" "<<phthis<<" "<<phnext<<"  dphi="<<dphi<<endl; 
      //cout<<"m,m',m'':"<<mthis<<" "<<mlast<<" "<<(-phlast*mnext+phnext*mlast + phthis*(mnext-mlast))/dphi<<endl;
      //dphi=phthis-phlast;
      //cout<<i<<" m ,Dm', Dm "<<mthis<<" "<<mthis-mlast<<" "<<Dm<<" ?>? "<<tolcut/dphi<<endl;
      //cout<<"  dphi="<< dphi<< endl;
      int i0=((i-1)%nphi+nphi)%nphi;
      //cout<<"i="<<i<<" i0="<<i0<<"  norefine="<<norefine<<endl;
      //if((not norefine) and fabs(Dm*dphi) > tolcut and dphi < 2*dphimin )cout<<"dphimin"<<endl;
      if((not norefine) and fabs(Dm) > tolcut and dphi >= 2*dphimin ){//fail test
	//cout<<"refine circle"<<endl;
	double dphidj=(phthis-phlast)/(nrefine+1.0);
	for(int j=-nrefine;j<0;j++){
	  double phi=phis[i0]+j*dphidj;
	  if(phi<0)phi+=twopi;
	  if(phi>=twopi)phi-=twopi;
	  newphis.push_back(phi);
	  newmags.push_back(0);
	}
	newphis.push_back(phis[i0]);
	newmags.push_back(mthis);
	dphidj=(phnext-phthis)/(nrefine+1.0);
	for(int j=1;j<=nrefine;j++){
	  double phi=phis[i0]+j*dphidj;
	  if(phi<0)phi+=twopi;
	  if(phi>=twopi)phi-=twopi;
	  newphis.push_back(phi);
	  newmags.push_back(0);
	}
	//cout<<i<<" "<<phlast<<" "<<phthis<<" "<<phnext<<"  dphi="<<dphi<<endl; 
	//cout<<"-->";
	//for(auto ph : newphis)cout<<ph<<" ";
	//cout<<endl<<"   ";
	//for(auto a : newmags)cout<<a<<" ";
	//cout<<endl;
	norefine=true; //Because we refine forward and backward of the current point.  It doesn't make sense to refine two points consecutively
	if(i==0)norefinelast=true;//If refining the first, then can't refine the last;
      } else {
	//cout<<"phthis/phi0:"<<phthis<<"/"<<phis[i0]<<endl;
	//cout<<"mthis/m0:"<<mthis<<"/"<<mags[i0]<<endl;
	newphis.push_back(phis[i0]);
	newmags.push_back(mthis);
	norefine=false;
      }
      mlast=mthis;
      mthis=mnext;
      phlast=phthis;
      phthis=phnext;
      //cout<<endl;
    }
    phis=newphis;
    mags=newmags;
    nphi=phis.size();
    //cout<<"Phis ";
    //for(auto ph : phis)cout<<ph<<" ";
    //cout<<endl<<"Mags ";
    //for(auto ph : mags)cout<<ph<<" ";
    //cout<<endl;
    //cout<<"circle: nphi="<<nphi<<"\n"<<endl;
  }
  //Now compute the mean mag
  double msum=0;
  for(int i=0;i<nphi;i++){
    int il=((i-1) % nphi + nphi) % nphi;
    double dphi=phis[i]-phis[il];
    if(dphi<0)dphi+=twopi;
    if(dphi>twopi)dphi-=twopi;
    msum+=(mags[i]+mags[il])*dphi/2.0;
  }
  //cout<<" r="<<radius<<" nphi="<<nphi<<" ctol="<<ctol<<endl;;
  magnification=msum/twopi;
  //cout<<"nphi="<<nphi<<endl;
  //cout<<"circle r="<<radius<<" mag="<<magnification<<"-----------\n"<<endl;
  return nphi;
}

///Brute force integration finite size magnification
///
///The idea here is that we perform a 2D integral over the image plane
int GLens::brute_force_area_mag(const Point &p, const double radius, double &magnification){
  const double tol=finite_source_tol;
  const double drhomin=tol*radius;
  const double tolcutmin=1e-9;
  const double maxctol=tol/30; //extremely slow below 1e5
  //const double maxctol=1e-4;  //compromise
  //const double maxctol=1e-3;
  //int nrho=4*(2+1/sqrt(tol));
  //const double maxctol=1e-3; //slightly faster
  int nrho=(20+1/sqrt(tol));
  //int nrho=3;
  int nrefine=2;
  const double maxmag=1e15;//for debug only

  //set up the initial grid in rho
  vector<double> mags(nrho,0);
  vector<double> rhos(nrho);
  for( int i=0;i<nrho;i++)rhos[i]=i*1.0/(nrho-1);

  //cout<<"\n********************\narea  radius="<<radius<<"  p="<<p.x<<" "<<p.y<<endl;
  double oldnrho=0;
  //main loop of refinement passes over the set of radii
  magnification=1;
  int neval=0;
  int met_tol_count=0;
  int icyc=0;
  while(oldnrho<nrho){
    //cout<<"p = "<<p.x<<" "<<p.y<<" "<<icyc<<endl;
    icyc++;
    oldnrho=nrho;
    //double margin=5+pow(nrho,.5);
    //double margin=5+pow(nrho,.5)/10;
    double margin=5+pow(nrho,.5)/3;
    //double margin=5+pow(nrho,.6667);
    //double cmargin=5+pow(nrho,.5);
    double cmargin=20;
    //double margin=sqrt(5+sqrt(nrho));
    //double tolcut=tol/margin;
    double tolcut=tol/margin;
    if(tolcut<tolcutmin)tolcut=tolcutmin;
    //cout<<"tolcut="<<tolcut<<endl;
    //compute necessary circle magnifications
    double ctolmin=1,ctolmax=0,ctolsum=0;//tuning diagnostics
    double tolcutmin=1,tolcutmax=0,tolcutsum=0;//tuning diagnostics
    double damin=1,damax=0,dasum=0;//tuning diagnostics
    int ctolcount=0;
    for( int i=0;i<nrho;i++){
      //cout<<"nrho,oldnrho,i: "<<nrho<<" "<<oldnrho<<" "<<i<<endl;
      if(mags[i]==0){
	double mag=0;
	if(i==0)mag=this->mag(invmap(p));
	else{
	  //the weighting for the tolerance on the circle mag estimate is based on the
	  //weighting farther below through which the circle mag estimate impact the area mag tolerance
	  double da=rhos[i]*rhos[i]-rhos[i-1]*rhos[i-1];
	  if(i<nrho-1)da=rhos[i+1]*rhos[i+1]-rhos[i-1]*rhos[i-1];
	  //double ctol=tolcut/da/sqrt(margin);//add a safety margin.
	  //double ctol=tolcut/da/cmargin;//add a safety margin.
	  double ctol=maxctol;
	  if(ctol<tolcutmin)ctol=tolcutmin;
	  if(ctol>maxctol)ctol=maxctol;
	  int evals=brute_force_circle_mag(p, rhos[i]*radius, ctol, mag);
	  neval+=evals;
	  if(0){//diagnostics
	    ctolcount++;
	    ctolsum+=log10(ctol);
	    if(ctol<ctolmin)ctolmin=ctol;
	    if(ctol>ctolmax)ctolmax=ctol;
	    cout<<rhos[i]<<" "<<mag<<" -- "<<evals<<" "<<log10(ctolmin)<<" < "<<log10(ctol)<<" , "<<ctolsum/ctolcount<<" < "<<log10(ctolmax)<<" < "<<log10(maxctol)<<endl;
	    tolcutsum+=log10(tolcut);
	    if(tolcut<tolcutmin)tolcutmin=tolcut;
	    if(tolcut>tolcutmax)tolcutmax=tolcut;
	    dasum+=log10(da);
	    if(da<damin)damin=da;
	    if(da>damax)damax=da;
	  }
	}
	//cout<<"i="<<i<<" mag="<<mag<<" p="<<p.x<<" "<<p.y<<" ->"<<this->mag(p)<<endl;
	double intens=1.0;//can make this a function of r for radial intensity profile
	mags[i]=mag*intens;
      }
    }
    if(0){
      int iprec = cout.precision();cout.precision(5);
      cout<<"log ctol: "<<log10(ctolmin)<<" < "<<ctolsum/ctolcount<<" < "<<log10(ctolmax);
      cout<<" log tolcut: "<<log10(tolcutmin)<<" < "<<tolcutsum/ctolcount<<" < "<<log10(tolcutmax)<<"  ";
      cout<<" log da: "<<log10(damin)<<" < "<<dasum/ctolcount<<" < "<<log10(damax)<<endl;
      cout.precision(iprec);
    }
    //cout<<endl;
    if(1){
      //Now compute the mean mag (for debugging)
      double msum=0;
      for(int i=1;i<nrho;i++){
	double da=(rhos[i]-rhos[i-1])*(rhos[i]+rhos[i-1]);//note area fraction adds identically to 1
	msum+=da*(mags[i]+mags[i-1])/2;
      }
      //cout<<"msum , msum-last :"<<msum<<" "<<msum-magnification<<endl;
      if(fabs(msum-magnification)<tol){
	if(met_tol_count>1)break;
	met_tol_count++;
      } else met_tol_count=0;
      magnification=msum;
    }
    
    //now test and construct the new grid (inelegant, but safe...)
    vector<double> newmags;
    vector<double> newrhos;
    double mll=0,ml=0,rl=0,rll=0;
    for(int i=0;i<nrho;i++){
      //cout<<"mags["<<i<<"]="<<mags[i]<<endl;
      double m0=mags[i];
      double r0=rhos[i];
      if(i>0){
	double dM;
	double r02=r0*r0,rl2=rl*rl,rll2=rll*rll;
	if(i==1){//err estimated from linear-m cubic-drho
	  //double drho=r0-rl;
	  //dM=(m0-ml)*drho*drho/6;...*(r0-rl)
	  dM=(m0-ml)*(r02-rl2);
	} else {//err est from quad-m cubic-rho
	  //cout<<"rll-r0: "<<rll<<" "<<rl<<" "<<r0<<endl;
	  //cout<<"mll-m0: "<<mll<<" "<<ml<<" "<<m0<<endl;
	  //dM=((m0-ml)/(r0-rl)-(ml-mll)/(rl-rll))*r0*r0*r0;...*(r0-rl)
	  dM=((m0-ml)*(r02-rl2)-(ml-mll)*(rl2-rll2))/2.0; //This is the change in net M from inserting point at rl between rll and r0
	  dM*=(r0-rl)/(r0-rll);//we are considering refining between rl and r0 not rll and r0, so reweight.
	}
	//cout<<i<<" DM = "<<dM<<endl; 
	//test magnitude variation tolerance
	//if(fabs(dM)*(r0-rl)*r0 > tolcut and (r0-rl) >= drhomin ){//fail test -> refine
	if(fabs(dM) > tolcut and (r0-rl) >= drhomin ){//fail test -> refine
	  //cout<<"refine disk"<<endl;
	  for(int j=0;j<nrefine;j++){
	    //cout<<"j="<<j<<"<"<<nrefine<<endl;
	    newrhos.push_back(rhos[i-1]+(rhos[i]-rhos[i-1])*(j+1.0)/(nrefine+1.0));
	    //cout<<"new rho="<<rhos[rhos.size()-1]<<endl;
	    newmags.push_back(0);
	  }
	}
      }
      rll=rl;mll=ml;
      rl=r0;ml=m0;	     
      newrhos.push_back(rhos[i]);
      newmags.push_back(mags[i]);
    }
    rhos=newrhos;
    mags=newmags;
    nrho=rhos.size();
    //cout<<"nrho="<<nrho<<endl;
  }
  //Now compute the mean mag
  double msum=0;
  for(int i=1;i<nrho;i++){
    double da=(rhos[i]-rhos[i-1])*(rhos[i]+rhos[i-1]);//note area fraction adds identically to 1
    msum+=da*(mags[i]+mags[i-1])/2;
  }
  magnification=msum;
  // cout<<"nrho="<<nrho<<"\n"<<endl;
  cout<<"nrho="<<nrho<<" neval="<<neval<<endl;
  cout<<"area  mag="<<magnification<<endl;
  //cout<<"******"<<endl;
  return neval;
}
  
///Brute mapping finite size magnification
///
///This method is intended only to be guarantedly correct, not fast. The idea is to map the relevant area of the image plane in a field of sample points, then to map these points back to the source plane to compute the magnification. This could also be a moderately efficient way to implement arbitrary source face intensity profiles.
int GLens::brute_force_map_mag(const Point &p, const double radius, double &magnification){

  //control parameters 
  double tol=finite_source_tol;
  double h=radius*sqrt(tol/M_PI);//assuming linear convergence in pixel area=h^2 
  int Npoly=1/sqrt(tol);
  const double box_expansion_fac=10;
  
  Point b=p;
  double Amag=0;
  double variance;
  const double radius2=radius*radius;
  vector< vector <Point> > contours;
  
  //First we call image area mag, needed to provide a location for the contours
  image_area_mag(b, radius, Npoly, Amag, variance, NULL, &contours);  
  
  //Make boxes around the computed image contours
  struct box {int xmin,xmax,ymin,ymax;};
  vector<box> boxes;
  vector<Point> ths=this->invmap(p);int ic=0;//diagn.
  for(auto curve:contours){
    //Find bounding box
    double xmin=1e100,ymin=1e100,xmax=-1e100,ymax=-1e100;
    for(auto p:curve){
      if(p.x<xmin)xmin=p.x;
      if(p.y<ymin)ymin=p.y;
      if(p.x>xmax)xmax=p.x;
      if(p.y>ymax)ymax=p.y;
    }
    //cout<<"countour box:"<<xmin<<" "<<xmax<<" "<<ymin<<" "<<ymax<<endl;
    //cout<<"point images:      "<<ths[ic].x<<"        "<<ths[ic].y<<endl;
    //ic++;
    //Expand boxes slightly, snapped to grid of resolution h    
    int ixc=int(((xmin+xmax)/2-p.x)/h+0.5);
    int iyc=int(((ymin+ymax)/2-p.y)/h+0.5);
    int nx=int(box_expansion_fac*(xmax-xmin)/h/2);//one-sided count n
    int ny=int(box_expansion_fac*(ymax-ymin)/h/2);//total points = 2*n+1
    //cout<<"box:"<<ixc-nx<<" "<<ixc+nx<<" "<<iyc-ny<<" "<<iyc+ny<<endl;
    boxes.push_back((box){ixc-nx,ixc+nx,iyc-ny,iyc+ny});
  }

  //Sum up weighted image area in grid units
  double sum=0;
  for(int k=0;k<boxes.size();k++){
    auto box=boxes[k];
    //cout<<"box "<<k<<endl;
    for(int i=box.xmin;i<=box.xmax;i++){
      double x=p.x+h*i;
      //check which boxes may already have counted this row
      int jcountedmin=box.ymax+1,jcountedmax=box.ymin-1;
      for(int kk=0;kk<k;kk++)
	if( i>=boxes[kk].xmin and i<=boxes[kk].xmax
	    and box.ymax>=boxes[kk].ymin and box.ymin<=boxes[kk].ymax){//overlaps in x
	  if(boxes[kk].ymax>=box.ymin and boxes[kk].ymin<jcountedmin)jcountedmin=boxes[kk].ymin;
	  if(boxes[kk].ymin<=box.ymax and boxes[kk].ymax>jcountedmax)jcountedmax=boxes[kk].ymax;
	  //cout<<"overlap range ["<<jcountedmin<<","<<jcountedmax<<"]"<<endl;
	}
      for(int j=box.ymin;j<=box.ymax;j++){
	//check if this cell is already counted
	if(j>=jcountedmin and j<=jcountedmax)continue;
	//map back to source plane to compute pixel intensity (uniform disk implemented)
	Point pgrid(x,p.y+h*j);
	Point psource=this->map(pgrid);
	//cout<<h<<": pgrid("<<pgrid.x<<","<<pgrid.y<<")"<<" psource("<<psource.x<<","<<psource.y<<")"<<" p("<<p.x<<","<<p.y<<")"<<endl;
	//compute squared distance to source center
	double dx=psource.x-p.x, dy=psource.y-p.y;
	double d2=dx*dx+dy*dy;
	//cout<<"  "<<j<<" "<<d2<<" "<<radius2<<endl;

	double edgebuf=0;
	if(0){//handle edges at higher orger
	  //when is (d-radius)^2<2 -> |1-d2/radius2|/(1+d/radius)<sqrt2
	  //                       ~> |1-d2/radius2|<2*sqrt(2)
	  //                       ~> -2*sqrt2 < d2/radius2-1 <2sqrt2
	  edgebuf=1.0-2.0*sqrt(2.0);
	}
	//linear inclusion in disk computation
	double area=0;
	if(d2/radius2<1-edgebuf)area=1.0;//inside source disk
	else if(fabs(d2/radius2-1)<edgebuf){//handle fractional part of edge cell to O(h^4)
	  if(dy==0){
	    area=(radius-fabs(dx))/h+0.5;
	  }else{//this has been though through only for d2>rho2 ... needs work.
	    double z=(d2-sqrt(d2)*radius)/h;
	    double yp=(-dx/2 + z)/dy, ym=(-dx/2 - z)/dy;
	    double xp=0.5,xm=-0.5;
	    if(yp>0.5){//out top to right
	      yp=0.5;
	      xp=(-dy/2 + z)/dx;
	      if(xp>0.5)xp=0.5;
	    } else if(yp<-0.5){//out bottom to right
	      yp=-0.5;
	      xp=(-dy/2 + z)/dx;
	      if(xp>0.5)xp=0.5;
	    }
	    if(ym<-0.5){//in from bottom on left
	      ym=-0.5;
	      xm=(-dy/2 - z)/dx;
	      if(xm<-0.5)xm=-0.5;
	    }else if(ym>0.5){//in from top on left
	      ym=0.5;
	      xm=(-dy/2 + z)/dx;
	      if(xm<-0.5)xm=-0.5;
	    }
	    area=0.5*(1-(0.5-yp)*(0.5-xp)-(ym+0.5)*(xm+0.5));
	  }
	}
	sum+=area;
      }
    }
  }
  magnification=sum*h*h/radius2/M_PI;
  double dmag=magnification-Amag;
  //cout<<"mag dmag:"<<magnification<<" "<<dmag<<endl;
}

vector<double> GLens::_compute_trajectory_dummy_dmag;//dummy argument

void GLens::finite_source_compute_trajectory (const Trajectory &traj, vector<double> &time_series, vector<vector<Point> > &thetas_series, vector<double>&mag_series, vector<double>&dmag_series, ostream *out){
  //Can optionally provide out stream to which to write image curves.

  //Controls
  const bool debug=false;
  //const int Npoly_max=finite_source_Npoly_max;        //Part of magnification-based estimate for polygon order.
  const int Npoly_max=finite_source_Npoly_max*(1+16*(source_radius<1?source_radius:1));  //Empirical hack based on limited example
  const double Npoly_Asat=2.0;   //Saturate at Npoly_max when image_area/pi = Npoly_Asat 
  bool dont_mix= false;
  bool do_laplacian = false;
  bool do_polygon = false;
  bool do_brute = false;

  const double rho2=source_radius*source_radius;
  //const double mtol=1e-5;
  //const double ftol=1e-2*source_radius;
  const double mtol=1e-6;
  const double ftol=1e-2*source_radius;
  const double mag_lcut=mtol/rho2;
  double mag_pcut=mtol/rho2/rho2;
  const double shear_cut=0.25/rho2/rho2;
  bool do_shear_test=false;
  //if(mag_pcut>1) mag_pcut=1.0;
  const double dmag_pcut=ftol;
  //decimation
  double decimate_dtmin=source_radius*finite_source_decimate_dtmin;

  //diagnostics
  bool diagnose=false;
  static int d_count=0,d_N=0,d_every=5000;
  static double d_time=0,d_rho,t_rho,t_N;
  double tstart=omp_get_wtime();
  int Nsum=0;

  int Npoly;

  if(finite_source_method<0)dont_mix=true;
  if(abs(finite_source_method)==1){
    do_polygon=true;
    do_laplacian=true;
  } else if(abs(finite_source_method)==2)do_laplacian=true;
  else if(abs(finite_source_method)==4)do_brute=true;
  
  if(debug){
    cout<<"source_radius="<<source_radius<<endl;
    cout<<"do_laplacian="<<do_laplacian<<endl;
    cout<<"do_polygon="<<do_polygon<<endl;
  }
  
  int Ngrid=traj.Nsamples();
  vector<double>full_time_series(Ngrid);
  vector<int>time_series_map;
  time_series.resize(0);
  double tlast=-INFINITY;
  for(int i=0; i<Ngrid;i++){
    double t=full_time_series[i]=traj.get_obs_time(i);
    if(i+2>=Ngrid or traj.get_obs_time(i+1)-tlast>=decimate_dtmin){
      time_series.push_back(t);
      time_series_map.push_back(i); 
      tlast=t;
    }
  }
  int Neval=time_series.size();
  //cout<<"Neval="<<Neval<<" < "<<Ngrid<<endl;
  thetas_series.resize(Neval);
  mag_series.resize(Neval);
  dmag_series.resize(Neval);
  
  for(int i=0; i<Neval;i++){
    double tgrid=time_series[i];
    Point b=get_obs_pos(traj,tgrid);
    if(debug)cout<<i<<" t="<<tgrid<<" b=("<<b.x<<","<<b.y<<")"<<endl;
    double Amag=0;
    Point CoM;
    double variance=0;

    
    //Here begins a decision tree of various possible finite source treatments

    //The first option is just to explicitly compute by brute force
    //This option works independently without mixing (A mix version could also be implemented below if desired
    //In this case,  we do not compute any variance or centroid information
    if(do_brute){
      //cout<<"t "<<tgrid<<endl;
      Npoly=brute_force_map_mag(b, source_radius, Amag);
      //Npoly=brute_force_area_mag(b, source_radius, Amag);
      //Pack up results
      mag_series[i]=Amag;
      dmag_series[i]=0;
      thetas_series[i]=vector<Point>(1,CoM-b);//TBD Except for the brute case, return the centriod. We haven't computed the centroid yet for brute.
      //cout<<Npoly<<" "<<tgrid<<" "<<Amag<<endl;
      Nsum+=Npoly;
      //cout<<i<<" mg0="<<mg0<<" Amag="<<Amag<<endl;
      continue; //finish the loop here. (Rest of this file is irrelevant in this case)
    }
    
    //At first we just compute the ordinary magnification and a leading-order finite source term
    //This is probably relatively fast enough that we can do it without worry about the additional cost
    vector<Point>thetas=invmap(b);
    int nk=thetas.size();
    double mg0 = mag(thetas);
    Nsum++;
    vector<double> mu0s(nk),mus(nk);
    for(int k=0;k<nk;k++)mu0s[k]=mag(thetas[k]);
    if(debug){
      cout<<"mg0="<<mg0<<endl;
      for(int k=0;k<nk;k++)cout<<"  mu0s["<<k<<"]="<<mu0s[k]<<endl;
    }
    //Now estimate the leading order finite source term:
    //This is a smaller calculation than the full Laplacian, keep in only up to 1/r^6 terms
    // dA*/A=1 + 4*norm(mu*dgamma)
    //Interestingly, for each of the images, the same term is dominant, and it is dA~O(1/r^6)
    //For the near-lens images, mu is small but dgamma is large, yielding the same order.
    //Analytically I get leading order for a binary as:
    //
    // A* = 1 + 2(1-q(1-q))r^(-4)( 1 + 4rho^2/r^2 )
    //
    //where all 3 images are included.
    bool shear_test=false;
    for(int k=0;k<nk;k++){
      Point th=thetas[k];
      vector<complex<double> > gammas;
      if(do_shear_test)gammas=compute_shear(th,2);
      else gammas=compute_shear(th,1);
      double dArel=0;
      if(mu0s[k]!=0)dArel=norm(mu0s[k]*gammas[1]*source_radius);
      double Ak=abs(mu0s[k])*(1+dArel);
      //cout<<"k="<<k<<gammas[1]<<" "<<dArel<<" "<<Ak<<endl;
      Amag+=Ak;
      CoM=CoM+th*Ak;
      //This test is based on an estimator for the max difference in mu^-1 near a point where mu^-1=0 only for the outer mu>1 image. 
      if(do_shear_test and mu0s[k]>1){
	cout<<"shear test: "<<shear_cut<<" < "<< norm(gammas[0]*gammas[2]) + norm(gammas[1]*gammas[1]) <<" mu="<<mu0s[k]<<" "<<mg0<<endl;
	shear_test = ( shear_cut < norm(gammas[0]*gammas[2]) + norm(gammas[1]*gammas[1]) );
      }
    }
    CoM=CoM*(1.0/Amag);
    if(debug)cout<<"Amag[lo]="<<Amag<<endl;
    if(debug)cout<<"CoM[lo]="<<CoM.x<<" "<<CoM.y<<endl;
    
    //Eventually we will want to dynamically select efficient methods for different regions.
    //For now we just have fixed choice of analytic or polygon methods
    if(debug)cout<<  Amag -1 <<" > "<<mag_lcut<<" ? "<<( Amag - 1 > mag_lcut)<<endl;
    bool do_laplacian_test= do_laplacian and ( Amag - 1 > mag_lcut or dont_mix);
    if(debug)cout<<" do_laplacian_test="<<do_laplacian_test<<endl;
    if(do_laplacian_test){
      if(debug)cout<<"doing laplacian"<<endl;
      //This method builds on PejchaEA2007? method
      // Amag = \sum_k I[k]/I[0] Lap^k[mu] / (2^k k!)^2
      // I[k] = \int_0^1 r^(2k+1) B(r) dr
      //Where B(rho/rho*) is the surface brightness at radius rho, for star-disk of radius rho*.
      //For constant surface brightness, I[k]/I[0]=1/(k+1)
      //Here we just keep the first nonleading term
      Amag=0;
      CoM=Point(0,0);
      for(int k=0;k<nk;k++){
	if(mu0s[k]==0)continue;
	Point th=thetas[k];
	double Lmu=Laplacian_mu(th);
	double dArel=Lmu*source_radius*source_radius/4.0/mu0s[k];
	dArel/=2.0;  //This is a total experimental HACK playing around...
	double Ak=abs(mu0s[k])*(1+dArel);
	if(debug)cout<<k<<" "<<th.x<<" "<<th.y<<" L="<<Lmu<<" mu="<<mu0s[k]<<" Ak="<<endl;
	Amag+=Ak;
	CoM=CoM+th*Ak;
      }
      CoM=CoM*(1.0/Amag);
    }

    
    if(debug)cout<<  Amag -1 <<" > "<<mag_pcut<<" ? or "<<  abs(Amag/mg0 - 1)  <<" > "<<dmag_pcut<<" ?  shear_test="<<shear_test<<endl;
    bool do_polygon_test= do_polygon  and ( shear_test or Amag - 1 > mag_pcut or abs(Amag/mg0 - 1) > dmag_pcut or dont_mix);
    if(debug)cout<<" do_polygon_test="<<do_polygon_test<<endl;
    if(do_polygon_test){
      if(debug)cout<<"doing polygon"<<endl;
      //This section computes the polygon order to apply
      //There are several possibilities in principle:
      //  -Use adaptive stepping in the polygon computation itself (maybe best long term)
      //  -Use an estimate from an analytic estimate of size of finite source effect
      //  -Use an estimate based on point-source magnification [implememted here]
      //  -Fixed (probably way too slow).
      //
      //Magnification-based Npoly: Based on the idea that the mean side length of poly is fixed
      //  -scales with sqrt(mg)
      //  -min of 4          as  Area/pi -> rho^2
      //  -max of Npoly_max  as  Area/pi >= Npoly_Asat
      //  -always even (to preserve time symmetry)
      const double N2scale=Npoly_max*Npoly_max/16.0-1.0;//4*sqrt(N2scale+1)=Npoly_max
      double extra_area = (mg0-1)/Npoly_Asat;
      if(extra_area>1.0)extra_area=1.0;       //extra_area ranges from 0 to 1
      //int Npoly = 2 * (int)(2*sqrt(1.0 + extra_area*N2scale));
      Npoly = 2 * (int)(2*sqrt(1.0 + extra_area*extra_area*N2scale));
      Npoly*=30;
      if(debug)cout<<" Npoly="<<Npoly<<" < "<<Npoly_max<<" N2scale="<<N2scale<<" extra_area="<<extra_area<<" Npoly="<<Npoly<<endl;
      //results go in Amag and CoM
      //cout<<" mu_i={ ";for(auto mui : mu0s)cout<<mui<<" ";cout<<"}"<<endl;
      //cout<<"source_radius="<<source_radius<<endl;
      image_area_mag(b, source_radius, Npoly, Amag, variance, out);  
      //if(Npoly>Npoly_max*5)cout<<"Npoly="<<Npoly<<endl;
      CoM=b;
      Nsum+=Npoly;
    }
    //Sanity check
    if(Amag<1){
      //if(1-Amag>1e-1)cout<<"impossible total magnification = "<<Amag<<" (polygon="<<do_polygon_test<<"), setting to mg0="<<mg0<<endl;
      Amag=mg0;
    }

    if(not do_polygon)cout<<" didn't do polygon"<<endl;
    //Pack it back up
    //cout<<"COM="<<CoM.x<<" "<<CoM.y<<endl;
    thetas_series[i]=vector<Point>(1,CoM-b); //Note we return lenght-1 vector with the overall image centroid offset.
    mag_series[i]=Amag;
    dmag_series[i]=sqrt(variance)*source_var;
    //cout<<i<<" mg0="<<mg0<<" Amag="<<Amag<<endl;

    //cout<<Npoly<<" "<<tgrid<<" "<<Amag<<endl;
    
    if(debug and do_polygon_test){
#pragma omp critical    
      cout<<"Npoly="<<Npoly<<" Amag="<<Amag<<" var="<<variance<<endl;
    }
  }

  if(debug){
#pragma omp critical
    cout<<"Nsum="<<Nsum<<endl;
  }
    
  if(Neval<Ngrid){//Need to reconstitute full grid by interpolation.
    //cout<<"Neval<Ngrid: "<<Neval<<"<"<<Ngrid<<endl;
    //cout<<"len(thetas)="<<thetas_series.size()<<endl;
    int imap=0;
    double t0,t1,m0,m1,v0,v1;
    Point c0,c1;
    vector<double>full_mag_series(Ngrid);
    vector<double>full_dmag_series(Ngrid);
    vector<vector<Point> > full_thetas_series(Ngrid);
    for(int i=0; i<Ngrid;i++){
      if(imap==0 or ( imap+1<time_series_map.size() and time_series_map[imap]<=i ) ){
	//cout<<"imap="<<imap<<endl;
	t0=time_series[imap];
	t1=time_series[imap+1];
	m0=mag_series[imap];
	m1=mag_series[imap+1];
	v0=dmag_series[imap];
	v1=dmag_series[imap+1];
	c0=thetas_series[imap][0];
	c1=thetas_series[imap+1][0];
	imap++;
      }
      double t=full_time_series[i];
      //fill values by linear interpolation, minimizing issues when edges are not smooth
      full_mag_series[i]=(m1*(t-t0)+m0*(t1-t))/(t1-t0);
      full_dmag_series[i]=(v1*(t-t0)+v0*(t1-t))/(t1-t0);
      vector<Point> thetas;
      thetas.push_back((c1*(t-t0)+c0*(t1-t))*(1.0/(t1-t0)));
      full_thetas_series[i]=thetas;
      /*
      cout<<"i="<<i<<endl;
      cout<<"  t: "<<t0<<", "<<t1<<" -> "<<full_time_series[i]<<endl;
      cout<<"  m: "<<m0<<", "<<m1<<" -> "<<full_mag_series[i]<<endl;
      cout<<"  v: "<<v0<<", "<<v1<<" -> "<<full_dmag_series[i]<<endl;
      cout<<" cx: "<<c0.x<<", "<<c1.x<<" -> "<<full_thetas_series[i][0].x<<endl;
      cout<<" cy: "<<c0.y<<", "<<c1.y<<" -> "<<full_thetas_series[i][0].y<<endl;
      */
    }
    /*
    cout<<"sizes time:"<<time_series.size()<<" -> "<<full_time_series.size()<<endl;
    cout<<"sizes  mag:"<<mag_series.size()<<" -> "<<full_mag_series.size()<<endl;
    cout<<"sizes dmag:"<<dmag_series.size()<<" -> "<<full_dmag_series.size()<<endl;
    cout<<"sizes thts:"<<thetas_series.size()<<" -> "<<full_thetas_series.size()<<endl;
    */

    //now copy back:
    mag_series=full_mag_series;
    dmag_series=full_dmag_series;
    thetas_series=full_thetas_series;
    time_series=full_time_series;
  }
  /*
  cout<<"sizes time:"<<time_series.size()<<endl;
  cout<<"sizes  mag:"<<mag_series.size()<<endl;
  cout<<"sizes dmag:"<<dmag_series.size()<<endl;
  cout<<"sizes thts:"<<thetas_series.size()<<endl;
  */

  if(diagnose){
    double dt=omp_get_wtime()-tstart;
#pragma omp critical
    {
      //if(source_radius>1)cout<<"rho="<<source_radius<<"!"<<endl;
      d_count++;
      d_N+=Nsum;
      d_rho+=source_radius;
      d_time+=dt;
      t_N+=Nsum*dt;
      t_rho+=source_radius*dt;
      if(d_count%d_every==0){
	cout<<"d_count="<<d_count<<" Nsum, rho, dt="<<Nsum<<" "<<source_radius<<" "<<dt<<endl;
	cout<<" <Nsum>, <rho>, <dt>="<<d_N/(double)d_count<<" "<<d_rho/d_count<<" "<<d_time/d_count<<endl;
	cout<<" <Nsum>_t, <rho>_t "<<t_N/(double)d_time<<" "<<t_rho/d_time<<endl;
      }
    }
  }
};
    
//Use GSL routine to integrate 
void GLens::compute_trajectory (const Trajectory &traj, vector<double> &time_series, vector<vector<Point> > &thetas_series, vector<int> &index_series,vector<double>&mag_series,vector<double> &dmag, bool integrate)
{
  // Given a trajectory through the observer plane, and a list of observation times, integrate the Jacobian to yield the corresponding trajectory in the lens plane.
  //
  //Arguments:
  //
  //Trajectory traj       -provides information about the trajectory of early thorough the observer plane.
  //traj->times  -provides a list of observation times to be included in the sample set 
  //vector<double> theta  -yields the resulting thetas for all sample times
  //vector<int> grid_idxs -has lenght of "times" and yield index values for the location of the corresponding results within the larger "theta" vector.
  //bool integrate (false for direct polynomial evaluation rather than integration. 
  //  if use_integrate is set then the value it overrides integrate 
  //
  //cout<<"lens="<<print_info()<<endl;
  //cout<<"compute_trajectory for traj="<<traj.print_info()<<endl;

  if(do_finite_source&&source_radius>0){//For finite-sources, we use a different approach
    //ostringstream oss;oss<<"curves_"<<source_radius<<".dat";
    //ofstream out(oss.str());
    //finite_source_compute_trajectory( traj, time_series, thetas_series, mag_series, dmag, &out);
    finite_source_compute_trajectory( traj, time_series, thetas_series, mag_series, dmag, finite_source_image_ofstream);
    int n=time_series.size();
    index_series.resize(n);
    for(int i=0;i<n;i++)index_series[i]=i;
    return;
  }
  
  //control parameters:
  double caustic_mag_poly_level = GL_int_mag_limit; //use direct polynomial eval near caustics.
  const double intTOL = GL_int_tol;  //control integration error tolerance
  const double test_result_tol = 1e-4;  //control integration error tolerance
  if(have_integrate)integrate=use_integrate;
  double prec=cout.precision();cout.precision(20);
  const double rWide_int_fac=100.0;

  
  //cout<<"glens::compTraj: int="<<integrate<<"\nthisLens="<<print_info()<<"\n traj="<<traj.print_info()<<endl;
  //cout<<"this="<<this<<endl;
  cout.precision(prec);

  ///clear the outputs
  time_series.clear();
  thetas_series.clear();
  index_series.clear();
  mag_series.clear();

  trajectory=&traj;//a convenience for passing to the integrator
  int NintSize=2*NimageMax;

  //Set up for GSL integration routines
  const gsl_odeiv2_step_type * stepType = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_step * step=nullptr;
  gsl_odeiv2_control * control=nullptr;
  gsl_odeiv2_evolve * evol=nullptr;
  gsl_odeiv2_system sys = {GLens::GSL_integration_func_vec, NULL, (size_t)NintSize, this};
  double h=1e-5;
  if(integrate){
    step = gsl_odeiv2_step_alloc (stepType, NintSize);
    control = gsl_odeiv2_control_y_new (intTOL, 0.0);
    //control = gsl_odeiv2_control_standard_new (intTOL, 0.0,1.0,10.0);//Also apply limits on the derivative
    evol = gsl_odeiv2_evolve_alloc (NintSize);
    gsl_odeiv2_control_init(control,intTOL,0.0,1,0);//Do I need this?
    gsl_odeiv2_evolve_reset(evol);//or this
  }

  //Main loop over observation times specified in the Trajectory object
  //initialization
  Point beta;
  vector<Point> thetas;
  bool evolving=false;
  double mg;
  have_saved_soln=false;

  int Ngrid=traj.Nsamples();
  double t_old=-1e100;
  for(int i=0; i<Ngrid;i++){    
    double tgrid=traj.get_obs_time(i);
    //entering main loop

    //We either compute the next step by evolution or by polynomial evaluation.
    //We do evolution is the integration flag is set *and* the last evaluated point had a
    //magnification under caustic_mag_poly_level.
    //If we evolving then the step-size may be smaller as driven by the ODE integrator's step-size 
    //control, we expect this to be small enough to avoide stumbling across the caustic, unaware.
    //After each step we check whether the mg level condition is satisfied, to consider whether to evolve.
    //We may also wish to implement a floor in the step-size near caustics for display purposes.

    if(evolving){
      double t=t_old;
      //The next loop steps (probably more finely) toward the next grid time point.      
      while(t<tgrid){
	//integrate each image
	if(debugint)cout<<"\n t="<<t<<" mg="<<mg<<endl;
	Ntheta=thetas.size();
	double theta[NintSize];
	for(int image=0;image<thetas.size();image++){
	  theta[2*image]=thetas[image].x;
	  theta[2*image+1]=thetas[image].y;
	}
	if(rWide_int_fac>0){
	  //Point b=traj.get_obs_pos(t);//TRAJ:Allow non-trivial transformation from observer-plane coords frame to lens frame coords
	  Point b=get_obs_pos(traj,t);
	  if(testWide(b,rWide_int_fac)){
	    //cout<<"evol: testWide==true"<<endl;
	    evolving=false;//switch to point-wise (WideBinary assuming rWide_int/rWide>=1)
	    break;
	  }
	}
	for(int k=2*thetas.size();k<NintSize;k++)theta[k]=0;
	int status = gsl_odeiv2_evolve_apply (evol, control, step, &sys, &t, tgrid, &h, theta);
	//Need some step-size control checking for near caustics?
	if (status != GSL_SUCCESS) { 
	  evolving=false;//switch to polynomial
	  break;
	}
	//record results
	for(int image=0;image<thetas.size();image++){
	  thetas[image]=Point(theta[2*image],theta[2*image+1]);	    
	}
	mg=mag(thetas);
	if(mg>=caustic_mag_poly_level)evolving=false;//switch to polynomial and don't record result
	else {
	  if(debugint)cout<<"mg="<<mg<<endl;
	  time_series.push_back(t);
	  thetas_series.push_back(thetas);
	  mag_series.push_back(mg);
	  if(!isfinite(mg)){
	    cout<<"integrate: mg is infinite!\n";
	    for(int image=0;image<thetas.size();image++){
	      cout<<theta[2*image]<<","<<theta[2*image+1]<<endl;
	    }	    
	  }
	}
      }
      if(!evolving){
	i--;//go back and try this step again with solving polynomial
	continue;
      }
    } else { //not evolving, solve polynomial
      beta=get_obs_pos(traj,tgrid);
      //cout<<i<<" t="<<tgrid<<" b=("<<beta.x<<","<<beta.y<<")"<<endl;
      //cout<<"Not evolving: beta=("<<beta.x<<","<<beta.y<<")"<<endl;
      thetas.clear();
      thetas=invmap(beta);
      //record results;
      mg=mag(thetas);
      time_series.push_back(tgrid);
      thetas_series.push_back(thetas);
      mag_series.push_back(mg);
      //cout<<"mg="<<mg<<", caustic_mag_poly_level="<<caustic_mag_poly_level<<endl;
      if(integrate&&mg<caustic_mag_poly_level){
	double r2=beta.x*beta.x+beta.y*beta.y;
	evolving=true;
	if(testWide(beta,rWide_int_fac))evolving=false;//Don't switch if in "wide" domain.
	//if(!evolving)cout<<"not evol: testWide==true"<<endl;
	if(!(thetas.size()==3||thetas.size()==5))evolving=false;//Don't switch to integrate if the number of images doesn't make sense
	if(evolving)gsl_odeiv2_evolve_reset(evol);
      };
      if(!isfinite(mg)){
	//We do some specific handling for GLensBinary here, could be cleaned up...
	GLensBinary* gb;
	bool squak=true;
	if(gb=dynamic_cast< GLensBinary* > (this)){
	  if(gb->get_q()<1e-14)squak=false;//we are going to fail with NAN, but not give a bunch of output, deal with it...
	  else cout<<"q,L="<<gb->get_q()<<","<<gb->get_L()<<endl;
	}
	if(squak){
	  cout<<"!integrate: mg is infinite! at beta="<<beta.x<<","<<beta.y<<endl;
	  for(int image=0;image<thetas.size();image++){
	    cout<<thetas[image].x<<","<<thetas[image].y<<" -> "<<mag(thetas[image])<<endl;
	  }
	  Trajectory::verbose=true;
	  beta=get_obs_pos(traj,tgrid);
	  alert();	 
	  Trajectory::verbose=false;
	}
      }
      if(debugint){
	cout<<"polynomial calc at t="<<tgrid<<":"<<endl;
	for(int image=0;image<thetas.size();image++)
	  cout<<"   theta["<<image<<"] = ("<<thetas[image].x<<","<<thetas[image].y<<")"<<endl;
      }
    }//end of polynomial step
    if(time_series.size()<1)cout<<"Time series empty i="<<i<<endl;
    t_old=tgrid;
    index_series.push_back(time_series.size()-1);
  }//end of main observation times loop

  if(integrate){
    gsl_odeiv2_evolve_free (evol);
    gsl_odeiv2_control_free (control);
    gsl_odeiv2_step_free (step);
  }

  if(test_result){
  //initialization
    for(int i=0; i<Ngrid;i++){
      double ttest=traj.get_obs_time(i);
      int ires=index_series[i];
      double tres=time_series[ires];
      Point beta=get_obs_pos(traj,ttest);
      vector<Point> thetas=invmap(beta);
      int nimages=thetas.size();
      double mgtest=mag(thetas);
      double mgres=mag_series[ires];
      double tminus,tplus,mgminus,mgplus;
      if(ires>0&&ires<time_series.size()-1){
	tminus=time_series[ires-1];
	tplus=time_series[ires+1];
	beta=get_obs_pos(traj,tminus);
	thetas=invmap(beta);
	mgminus=mag(thetas);
	beta=get_obs_pos(traj,tplus);
	thetas=invmap(beta);
	mgplus=mag(thetas);
      }

      if(abs(mgtest-mgres)/mgtest>test_result_tol)
#pragma omp critical
	{
	cout<<"\ncompute_trajectory: test failed. test/res:\nindex="<<i<<","<<ires<<"\ntime="<<ttest<<","<<tres<<"\nmag="<<mgtest<<","<<mgres<<" -> "<<abs(mgtest-mgres)<<"["<<nimages<<" images]"<<endl;
	cout<<"beta=("<<beta.x<<","<<beta.y<<")"<<endl;
	if(ires>0&&ires<time_series.size()-1){
	  cout<<"nearby times  :"<<tminus<<" < t < "<<tplus<<endl;
	  cout<<"   with mags  :"<<mgminus<<" < mg < "<<mgplus<<endl;
	}
	//if(i>0&&ires<index_series.size()-1)cout<<"nearby idx times:"<<time_series[index_series[i-1]]<<" < t < "<<time_series[index_series[i+1]]<<endl;
      }
    }
  }

}



/* Under development
//Use GSL routine to integrate polynomial roots
void GLens::compute_trajectory_integrate_roots (const Trajectory &traj, vector<double> &time_series, vector<vector<Point> > &thetas_series, vector<int> &index_series,vector<double>&mag_series,vector<double> &dmag)
{
  // Given a trajectory through the observer plane, and a list of observation times, integrate the Jacobian to yield the corresponding trajectory in the lens plane.
  //
  //Arguments:
  //
  //Trajectory traj       -provides information about the trajectory of early thorough the observer plane.
  //traj->times  -provides a list of observation times to be included in the sample set 
  //vector<double> theta  -yields the resulting thetas for all sample times
  //vector<int> grid_idxs -has lenght of "times" and yield index values for the location of the corresponding results within the larger "theta" vector.
  //bool integrate (false for direct polynomial evaluation rather than integration. 
  //  if use_integrate is set then the value it overrides integrate 
  //
  //cout<<"lens="<<print_info()<<endl;
  //cout<<"compute_trajectory for traj="<<traj.print_info()<<endl;

  //control parameters:
  double caustic_mag_poly_level = GL_int_mag_limit; //use direct polynomial eval near caustics.
  const double intTOL = GL_int_tol;  //control integration error tolerance
  const double test_result_tol = 1e-4;  //control integration error tolerance
  if(have_integrate)integrate=use_integrate;
  double prec=cout.precision();cout.precision(20);
  const double rWide_int_fac=100.0;

  cout.precision(prec);

  ///clear the outputs
  time_series.clear();
  thetas_series.clear();
  index_series.clear();
  mag_series.clear();

  trajectory=&traj;//a convenience for passing to the integrator
  int NintSize=2*NimageMax;

  //Set up for GSL integration routines
  const gsl_odeiv2_step_type * stepType = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_step * step=nullptr;
  gsl_odeiv2_control * control=nullptr;
  gsl_odeiv2_evolve * evol=nullptr;
  gsl_odeiv2_system sys = {GLens::GSL_integration_func_vec, NULL, (size_t)NintSize, this};
  double h=1e-5;
  step = gsl_odeiv2_step_alloc (stepType, NintSize);
  control = gsl_odeiv2_control_y_new (intTOL, 0.0);
  //control = gsl_odeiv2_control_standard_new (intTOL, 0.0,1.0,10.0);//Also apply limits on the derivative
  evol = gsl_odeiv2_evolve_alloc (NintSize);
  gsl_odeiv2_control_init(control,intTOL,0.0,1,0);//Do I need this?
  gsl_odeiv2_evolve_reset(evol);//or this


  //Main loop over observation times specified in the Trajectory object
  //initialization
  Point beta;
  vector<Point> thetas;
  bool evolving=false;
  double mg;
  have_saved_soln=false;

  int Ngrid=traj.Nsamples();
  double t_old=-1e100;
  for(int i=0; i<Ngrid;i++){    
    double tgrid=traj.get_obs_time(i);
    //entering main loop
    
    //test whether to evolve based on how close to "real" root pairs are...TBD
    //also test for wide binary (at target time) and test "have_roots"
    if(evolve){
      double t=t_old;
      //The next loop steps (probably more finely) toward the next grid time point.      
      while(t<tgrid){
	//integrate each root
	vector<Point>roots;
	double theta[NintSize];
	for(int iroot=0;iroot<roots.size();iroot++){
	  theta[2*iroot]=thetas[iroot].x;
	  theta[2*iroot+1]=thetas[iroot].y;
	}
	for(int k=2*thetas.size();k<NintSize;k++)theta[k]=0;
	int status = gsl_odeiv2_evolve_apply (evol, control, step, &sys, &t, tgrid, &h, theta);
	//Need some step-size control checking for near caustics?
	if (status != GSL_SUCCESS) { 
	  evolving=false;//switch to polynomial
	  break;
	}
	//record results
	for(int image=0;image<thetas.size();image++){
	  thetas[image]=Point(theta[2*image],theta[2*image+1]);	    
	}
	
	//Now test again if that there are not  near real root pairs
	//if it is s real root, then we can also check that magnification
	//has not changed sign.  Does that generalize to ghost roots??
	if(fail test) evolve=false;
	else{
	  time_series.push_back(t);
	  thetas_series.push_back(thetas);
	  mag_series.push_back(mg);//need to compute mg??
	  if(!isfinite(mg)){
	    cout<<"integrate: mg is infinite!\n";
	    for(int image=0;image<thetas.size();image++){
	      cout<<theta[2*image]<<","<<theta[2*image+1]<<endl;
	    }	    
	  }
	}
      }
      if(!evolving){
	i--;//go back and try this step again with solving polynomial
	continue;
      }
    } else { //not evolving, solve polynomial
      beta=get_obs_pos(traj,tgrid);
      //cout<<i<<" t="<<tgrid<<" b=("<<beta.x<<","<<beta.y<<")"<<endl;
      //cout<<"Not evolving: beta=("<<beta.x<<","<<beta.y<<")"<<endl;
      thetas.clear();
      thetas=invmap(beta);//Want to use wittmag with nocheck directly instead to get all root
      //record results?? First need to check roots tho...
      mg=mag(thetas);
      time_series.push_back(tgrid);
      thetas_series.push_back(thetas);
      mag_series.push_back(mg);
      //cout<<"mg="<<mg<<", caustic_mag_poly_level="<<caustic_mag_poly_level<<endl;
      if(integrate&&mg<caustic_mag_poly_level){
	double r2=beta.x*beta.x+beta.y*beta.y;
	evolving=true;
	if(testWide(beta,rWide_int_fac))evolving=false;//Don't switch if in "wide" domain.
	//if(!evolving)cout<<"not evol: testWide==true"<<endl;
	if(!(thetas.size()==3||thetas.size()==5))evolving=false;//Don't switch to integrate if the number of images doesn't make sense
	if(evolving)gsl_odeiv2_evolve_reset(evol);
      };
      if(!isfinite(mg)){
	//We do some specific handling for GLensBinary here, could be cleaned up...
	GLensBinary* gb;
	bool squak=true;
	if(gb=dynamic_cast< GLensBinary* > (this)){
	  if(gb->get_q()<1e-14)squak=false;//we are going to fail with NAN, but not give a bunch of output, deal with it...
	  else cout<<"q,L="<<gb->get_q()<<","<<gb->get_L()<<endl;
	}
	if(squak){
	  cout<<"!integrate: mg is infinite! at beta="<<beta.x<<","<<beta.y<<endl;
	  for(int image=0;image<thetas.size();image++){
	    cout<<thetas[image].x<<","<<thetas[image].y<<" -> "<<mag(thetas[image])<<endl;
	  }
	  Trajectory::verbose=true;
	  beta=get_obs_pos(traj,tgrid);
	  alert();	 
	  Trajectory::verbose=false;
	}
      }
      if(debugint){
	cout<<"polynomial calc at t="<<tgrid<<":"<<endl;
	for(int image=0;image<thetas.size();image++)
	  cout<<"   theta["<<image<<"] = ("<<thetas[image].x<<","<<thetas[image].y<<")"<<endl;
      }
    }//end of polynomial step
    if(time_series.size()<1)cout<<"Time series empty i="<<i<<endl;
    t_old=tgrid;
    index_series.push_back(time_series.size()-1);
  }//end of main observation times loop

  if(integrate){
    gsl_odeiv2_evolve_free (evol);
    gsl_odeiv2_control_free (control);
    gsl_odeiv2_step_free (step);
  }

  if(test_result){
  //initialization
    for(int i=0; i<Ngrid;i++){
      double ttest=traj.get_obs_time(i);
      int ires=index_series[i];
      double tres=time_series[ires];
      Point beta=get_obs_pos(traj,ttest);
      vector<Point> thetas=invmap(beta);
      int nimages=thetas.size();
      double mgtest=mag(thetas);
      double mgres=mag_series[ires];
      double tminus,tplus,mgminus,mgplus;
      if(ires>0&&ires<time_series.size()-1){
	tminus=time_series[ires-1];
	tplus=time_series[ires+1];
	beta=get_obs_pos(traj,tminus);
	thetas=invmap(beta);
	mgminus=mag(thetas);
	beta=get_obs_pos(traj,tplus);
	thetas=invmap(beta);
	mgplus=mag(thetas);
      }

      if(abs(mgtest-mgres)/mgtest>test_result_tol)
#pragma omp critical
	{
	cout<<"\ncompute_trajectory: test failed. test/res:\nindex="<<i<<","<<ires<<"\ntime="<<ttest<<","<<tres<<"\nmag="<<mgtest<<","<<mgres<<" -> "<<abs(mgtest-mgres)<<"["<<nimages<<" images]"<<endl;
	cout<<"beta=("<<beta.x<<","<<beta.y<<")"<<endl;
	if(ires>0&&ires<time_series.size()-1){
	  cout<<"nearby times  :"<<tminus<<" < t < "<<tplus<<endl;
	  cout<<"   with mags  :"<<mgminus<<" < mg < "<<mgplus<<endl;
	}
	//if(i>0&&ires<index_series.size()-1)cout<<"nearby idx times:"<<time_series[index_series[i-1]]<<" < t < "<<time_series[index_series[i+1]]<<endl;
      }
    }
  }

}
*/

//This version seems no longer used 04.02.2016..
int GLens::GSL_integration_func (double t, const double theta[], double thetadot[], void *instance){
  //We make the following cast static, thinking that that should realize faster integration
  //However, since we call functions which are virtual, this may/will give the wrong result if used
  //with a derived class that has overloaded those functions {get_obs_pos, get_obs_vel, invjac, map}
  GLens *thisobj = static_cast<GLens *>(instance);
  const Trajectory *traj=thisobj->trajectory;
  Point p(theta[0],theta[1]);
  double j00i,j10i,j01i,j11i,invJ = thisobj->invjac(p,j00i,j01i,j10i,j11i);
  Point beta0=thisobj->get_obs_pos(*traj,t);
  Point beta=thisobj->map(p);
  double dx=beta.x-beta0.x,dy=beta.y-beta0.y;
  Point betadot=thisobj->get_obs_vel(*traj,t);
  double adjbetadot[2];
  //double rscale=thisobj->estimate_scale(Point(theta[0],theta[1]));
  //if our image is near one of the lenses, then we need to tread carefully
  double rk=thisobj->kappa;//*rscale*0;
  adjbetadot[0] = betadot.x-rk*dx;
  adjbetadot[1] = betadot.y-rk*dy;
  thetadot[0] = j00i*adjbetadot[0]+j01i*adjbetadot[1];
  thetadot[1] = j10i*adjbetadot[0]+j11i*adjbetadot[1];

  if(!isfinite(invJ))cout<<"GLens::GSL_integration_func: invJ=inf, |beta|="<<sqrt(beta0.x*beta0.x+beta0.y*beta0.y)<<endl;

  return isfinite(invJ)?GSL_SUCCESS:3210123;
}


int GLens::GSL_integration_func_vec (double t, const double theta[], double thetadot[], void *instance){
  //We make the following cast static, thinking that that should realize faster integration
  //However, since we call functions which are virtual, this may/will give the wrong result if used
  //with a derived class that has overloaded those functions {get_obs_pos, get_obs_vel, invjac, map}
  GLens *thisobj = static_cast<GLens *>(instance);
  const Trajectory *traj=thisobj->trajectory;
  Point beta0=thisobj->get_obs_pos(*traj,t);
  Point betadot=thisobj->get_obs_vel(*traj,t);
  bool fail=false;
  for(int k=0;k<2*thisobj->NimageMax;k++)thetadot[k]=0;
  //cout<<"beta0=("<<beta0.x<<","<<beta0.y<<")"<<endl;
  for(int image =0; image<thisobj->Ntheta;image++){
    Point p(theta[image*2+0],theta[image*2+1]);
    double j00i,j10i,j01i,j11i,invJ = thisobj->invjac(p,j00i,j01i,j10i,j11i);
    Point beta=thisobj->map(p);
    double dx=beta.x-beta0.x,dy=beta.y-beta0.y;
    double adjbetadot[2];
    //double rscale=thisobj->estimate_scale(Point(theta[0],theta[1]));
    //if our image is near one of the lenses, then we need to tread carefully
    double rk=thisobj->kappa;//*rscale;
    adjbetadot[0] = betadot.x;
    adjbetadot[1] = betadot.y;
    if(isfinite(dx))adjbetadot[0] += -rk*dx;
    if(isfinite(dy))adjbetadot[1] += -rk*dy;
    if(isfinite(invJ)){
      thetadot[2*image]   = (j00i*adjbetadot[0]+j01i*adjbetadot[1]);
      thetadot[2*image+1] = (j10i*adjbetadot[0]+j11i*adjbetadot[1]);
    }
    if(!isfinite(invJ)){    
      //cout<<"GLens::GSL_integration_func_vec: FAIL"<<endl;
      //fail=true;
      //here we panic and rather than return NAN, we evolve as if the lensing effect is trivial...
      cout<<"GLens::GSL_integration_func_vec: invJ=inf, |beta|="<<sqrt(beta0.x*beta0.x+beta0.y*beta0.y)<<endl;
      thetadot[2*image]   = betadot.x;
      thetadot[2*image+1] = betadot.y;
    }
  }
  for(int k=2*thisobj->Ntheta;k<2*thisobj->NimageMax;k++)thetadot[k]=0;

  return !fail?GSL_SUCCESS:3210123;
}

void GLens::addOptions(Options &opt,const string &prefix){
  Optioned::addOptions(opt,prefix);
  //addTypeOptions(opt);
  opt.add(Option("GL_poly","Don't use integration method for lens magnification, use only the polynomial method."));
  //opt.add(Option("poly","Same as GL_poly for backward compatibility.  (Deprecated)"));
  opt.add(Option("GL_int_tol","Tolerance for GLens inversion integration. (1e-10)","1e-10"));
  opt.add(Option("GL_int_mag_limit","Magnitude where GLens inversion integration reverts to poly. (1.5)","1.5"));
  opt.add(Option("GL_int_kappa","Strength of driving term for GLens inversion. (0.1)","0.1"));
  opt.add(Option("GL_finite_source","Flag to turn on finite source fitting. Optional argument to provide method [leading,laplacian,polygon,(no arg default), uses fastest appropriate, up to specification or use eg 'strict_polygon']"));
  opt.add(Option("GL_finite_source_Npoly_max","Max number of sides in polygon source approximation.(10 default)","10"));
  opt.add(Option("GL_finite_source_var","Factor (roughly) for variance in surface brightness from uniformity.(0.01 default)","0.01"));
  opt.add(Option("GL_finite_source_log_rho_max","Set max uniform prior range for log_rho. (-100->gaussian prior default)","-100"));
  opt.add(Option("GL_finite_source_log_rho_min","Set min if uniform prior for log_rho. (-6.0 default)","-6"));
  opt.add(Option("GL_finite_source_refine_limit","Maximum refinement factor. (100.0 default)","100.0"));
  opt.add(Option("GL_finite_source_tol","Magnitude tolerance target. (1e-4 default)","1e-5"));
  opt.add(Option("GL_finite_source_decimate_dtmin","Interpolate time-steps closer than this fraction of source size. (default sqrt(GL_finite_source_tol))","-1"));
};

void GLens::setup(){
  set_integrate(!optSet("GL_poly")&&!optSet("poly"));
  *optValue("GL_int_tol")>>GL_int_tol;
  *optValue("GL_int_mag_limit")>>GL_int_mag_limit;
  *optValue("GL_int_kappa")>>kappa;
  double finite_source_log_rho_max;
  double finite_source_log_rho_min;
  if(optSet("GL_finite_source")){
    do_finite_source=true;
    *optValue("GL_finite_source_var")>>source_var;
    string method;
    *optValue("GL_finite_source")>>method;
    finite_source_method=0;
    if(method=="polygon"||method=="true")finite_source_method=1;
    else if(method=="laplacian")finite_source_method=2;
    else if(method=="leading")finite_source_method=3;
    else if(method=="strict_polygon")finite_source_method=-1;
    else if(method=="strict_laplacian")finite_source_method=-2;
    else if(method=="strict_brute")finite_source_method=-4;
    else{
      cout<<"GLens::setup: Finite source method '"<<method<<"' not recognized."<<endl;
      exit(1);
    }
    cout<<"finite source method = '"<<method<<"' -> "<<finite_source_method<<endl; 
    *optValue("GL_finite_source_Npoly_max")>>finite_source_Npoly_max;
    *optValue("GL_finite_source_log_rho_max")>>finite_source_log_rho_max;
    *optValue("GL_finite_source_log_rho_min")>>finite_source_log_rho_min;
    *optValue("GL_finite_source_refine_limit")>>finite_source_refine_limit;
    *optValue("GL_finite_source_decimate_dtmin")>>finite_source_decimate_dtmin;
    *optValue("GL_finite_source_tol")>>finite_source_tol;
    if(finite_source_decimate_dtmin<0)finite_source_decimate_dtmin=sqrt(finite_source_tol);
  }
  haveSetup();
  cout<<"GLens set up with:\n\tintegrate=";
  if(use_integrate)cout<<"true\n\tGL_int_tol="<<GL_int_tol<<"\n\tkappa="<<kappa<<endl;
  else cout<<"false"<<endl;
  if(do_finite_source){
    string names[] =                                      {"log_rho_star"};
    nativeSpace=stateSpace(1);
    nativeSpace.set_names(names);
    GLSpace=nativeSpace;
    const int uni=mixed_dist_product::uniform, gauss=mixed_dist_product::gaussian, pol=mixed_dist_product::polar; 
    valarray<double>    centers((initializer_list<double>){ -4.0  });
    valarray<double>     scales((initializer_list<double>){  1.0  });
    valarray<int>         types((initializer_list<int>)   { gauss });
    if(finite_source_log_rho_max>-100.0){//set uniform prior for log_rho
      centers[0]=(finite_source_log_rho_max+finite_source_log_rho_min)/2.0;
      scales[0]=(finite_source_log_rho_max-finite_source_log_rho_min)/2.0;
      types[0]= uni;
    }
    setPrior(new mixed_dist_product(&GLSpace,types,centers,scales));
  } else {
    nativeSpace=stateSpace(0);
    setPrior(new sampleable_probability_function(&nativeSpace));//dummy
  }
};

///Compute the Laplacian of the local image magnification explicitly
///
/// Using C-R (Wirtinger) calculus, the calculation is built on:
///
///    mu = ( 1 - gamma*gammac )^(-1)
/// gamma = \sum_i^N nu_i / (zc*zc)   =  dbetac/dz
///
/// Where postscript c indicates conjugate
///
/// Note that dgamma/dzc = 0.
///
/// In the manifestly flat source/observer beta plane, but transformed
/// to Wirtinger variables:
///
///   ds^2 = dz dzc
///  Lap f = 4 (d/dz)(d/dzc) f
///
/// Then, after diffeomorphism to lens plane [away from caustics]:
/// 
/// Lap f = mu \partial_a ( ginv^{ab}/mu \partial b f )
///    mu = srt(det[g]), as given above
///  ginv = mu^2 (    1    -gamma )
///              ( -gammac    1   )
/// after algebra:
///
/// Lap f = 2 mu [-gammac (d/dz)^2 + (2-1/mu)(d/dz)(d/dzc) -gamma (d/dzc)^2 + gamma dgammac/dz (d/dz) ] f
///
/// Lap mu = mu^3 ( -4 gammac^2( 3 mu^2 gammac (dgamma/dz)^2 + mu (d^2gamma/dz^2) )
///                 +2 ( 6 mu^2 - 6 mu + 1 )
///                 + CC 
///               )
/// The same approach can be applied again to get Lap^2 mu, with a somewhat more complicated result
/// involving up to d4gamma.
double GLens::Laplacian_mu(const Point &p)const{
  //Computuing the shear and derives is specialized to the lens-type, the rest is general...
  vector<complex<double> >gammas=compute_shear(p,2);
  complex<double>   gamma=gammas[0], gammac=conj(gamma), dgamma=gammas[1], d2gamma=gammas[2];
  double invmu = 1-norm(gamma);
  if(invmu==0)return 0;//Fail gracefully;  There is no appropriate result near caustics
  double mu=1/invmu, mu2=mu*mu;
  double term1 = real(-2.0*gammac*gammac*( 3.0*mu2*gammac*dgamma*dgamma + mu*d2gamma ));
  double term2 = real(norm(dgamma)*( 6*mu*(mu-1) + 1 ));
  return 4*mu*mu2*(term1+term2);
};
 
///Compute the complex lens shear, and some number of its derivatives
/// gamma = \sum_i^N nu_i / (zc*zc)   =  dbetac/dz
vector<complex<double> > GLens::compute_shear(const Point &p, int nder)const{
  double x=p.x,y=p.y;
  complex<double> z(x,y);
  vector<complex<double> >gammas;   
  complex<double>   gamma=1.0/z/z, gammac=conj(gamma);
  complex<double> dNgamma=gamma;
  gammas.push_back(gamma);
  for(int n=0;n<nder;n++){
    dNgamma *=-(n+2.0)/z;
    gammas.push_back(dNgamma);
  }
  return gammas;
};

//
// ******************************************************************
// GLensBinary routines *********************************************
// ******************************************************************
//

GLensBinary::GLensBinary(double q,double L,double phi0):q(q),L(L),phi0(phi0),sin_phi0(sin(phi0)),cos_phi0(cos(phi0)){
  typestring="GLens";
  option_name="BinaryLens";
  option_info="Fixed binary point-mass lens";
  NimageMax=5;
  NimageMin=3;
  nu=1/(1+q);
  rWide=5;
  do_remap_q=false;
  q_ref=0;
  idx_q=idx_L=idx_phi0=-1;
};

void GLensBinary::setup(){
  set_integrate(!optSet("GL_poly")&&!optSet("poly"));
  *optValue("GL_int_tol")>>GL_int_tol;
  *optValue("GL_int_mag_limit")>>GL_int_mag_limit;
  *optValue("GL_int_kappa")>>kappa;
  *optValue("GLB_rWide")>>rWide;
  haveSetup();
  //cout<<"GLens set up with:\n\tintegrate=";
  //if(use_integrate)cout<<"true\n\tGL_int_tol="<<GL_int_tol<<"\n\tkappa="<<kappa<<endl;
  //else cout<<"false"<<endl;
  double q0_val;
  *optValue("q0")>>q0_val;
  if(optSet("remap_q")){
    remap_q(q0_val);
  }
  GLens::setup();
  //set nativeSpace
  stateSpace space(3);
  string names[] =                                      {"logq","logL","phi0"};
  const int uni=mixed_dist_product::uniform, gauss=mixed_dist_product::gaussian, pol=mixed_dist_product::polar; 
  valarray<double>    centers((initializer_list<double>){   0.0,   0.0,  M_PI});
  valarray<double> halfwidths((initializer_list<double>){   log10(q0_val),   1.0,  M_PI});
  valarray<int>         types((initializer_list<int>)   {   uni, gauss,   uni});
  if(do_remap_q)names[3]="s(1+q)";
  space.set_bound(2,boundary(boundary::wrap,boundary::wrap,0,2*M_PI));//set 2-pi-wrapped space for phi0.
  space.set_names(names);  

  //cout<<"native::space=\n"<<nativeSpace.show()<<endl;
  //cout<<"native is:\n"<<nativePrior->show()<<endl;
  //cout<<"binary::space=\n"<<space.show()<<endl;
  GLBinarySpace=space;
  nativeSpace.attach(space);
  //cout<<"new native::space=\n"<<nativeSpace.show()<<endl;
  if(optSet("GLB_gauss_q"))types[0]=gauss;
  if(do_remap_q){
    double qq=2.0/(q_ref+1.0);
    double ds=0.5/(1.0+qq*qq); //ds=(1-s(q=1))/2
    centers[0]=1.0-ds;
    halfwidths[0]=ds;          //ie range=[s(q=1),s(q=inf)=1.0]
    types[0]=uni;
  }
  //setPrior(new mixed_dist_product(&nativeSpace,types,centers,halfwidths));
  binaryPrior=make_shared<mixed_dist_product>(&GLBinarySpace,types,centers,halfwidths);
  parentPrior=nativePrior;
  //cout<<"parent is:\n"<<parentPrior->show()<<endl;
  //cout<<"binary is:\n"<<binaryPrior->show()<<endl;
  //cout<<"parent space is:\n"<<parentPrior->get_space()->show()<<endl;
  //cout<<"binary space is:\n"<<binaryPrior->get_space()->show()<<endl;
  setPrior(new independent_dist_product(&nativeSpace,parentPrior.get(),binaryPrior.get()));//append to prior
  //cout<<"net prior is:\n"<<nativePrior->show()<<endl;
};

Point GLensBinary::map(const Point &p){
  ldouble x=p.x,y=p.y,x1=x-ldouble(L)/2.0L,x2=x+ldouble(L)/2.0L,r1sq=x1*x1+y*y,r2sq=x2*x2+y*y;
  ldouble c1=(1.0L-nu)/r1sq,c2=nu/r2sq;
  return Point(x-x1*c1-x2*c2,y-y*(c1+c2));
};

vector<Point> GLensBinary::invmap(const Point &p){
  const double rTest=1.1*rWide;
  double r2=p.x*p.x+p.y*p.y;
  if(testWide(p,1.0)){
    if(debug||inv_test_mode&&debugint){
      debug=true;
      cout<<"wide"<<endl;
    }
    vector<Point>thWB= invmapWideBinary(p);
    //if fails to converge (rare) revert to WittMao:
    if(thWB.size()==0){
      //cout<<"WideBinary failed to converge"<<endl;
      have_saved_soln=false;//Can't rely on save soln when jumping from WideBinary
      return invmapWittMao(p);
    }
    if(inv_test_mode&&r2<rTest*rTest){
      if(debug||debugint)cout<<"wide-test"<<endl;
      vector<Point>thWM= invmapWittMao(p);
      double mag_WB=mag(thWB),mag_WM=mag(thWM);
      if(abs(mag_WB-mag_WM)/mag_WM>1e-6){//the two methods don't agree.
	if(debug||debugint)cout<<"wide-test FAILED"<<endl;
	double TOL=LEADTOL*LEADTOL*(100000000+p.x*p.x+p.y*p.y+4.0/L/L);
	cout.precision(15);
	cout<<"\nInversion methods disagree: magWB="<<mag_WB<<" magWM="<<mag_WM<<" at ("<<p.x<<","<<p.y<<")"<<endl;
	cout<<"WittMaoCalc:"<<endl;
	debug=true;invmapWittMao(p);debug=false;
	cout<<"----"<<endl;
	for(int i=0;i<max(thWM.size(),thWB.size());i++){
	  if(i<thWB.size()){
	    Point th=thWB[i];
	    Point beta=map(th);
	    double dx=beta.x-p.x,dy=beta.y-p.y;
	    double dbx=th.x-p.x,dby=th.y-p.y;
	    cout<<"WB("<<i<<"): "<<th.x<<" "<<th.y<<" "<<mag(th)<<" : ("<<beta.x<<"-"<<p.x<<","<<beta.y<<"-"<<p.y<<")^2 "<<(dx*dx+dy*dy<TOL*(1+dbx*dbx+dby*dby)?" < ":"!< ")<<TOL*(1+dbx*dbx+dby*dby)<<endl;
	  } else {
	    cout<<"WB(---)"<<endl;
	  }
	  if(i<thWM.size()){
	    Point th=thWM[i];
	    Point beta=map(th);
	    double dx=beta.x-p.x,dy=beta.y-p.y;
	    double dbx=th.x-p.x,dby=th.y-p.y;
	    cout<<"WM("<<i<<"): "<<th.x<<" "<<th.y<<" "<<mag(th)<<" : ("<<beta.x<<"-"<<p.x<<","<<beta.y<<"-"<<p.y<<")^2 "<<(dx*dx+dy*dy<TOL*(1+dbx*dbx+dby*dby)?" < ":"!< ")<<TOL*(1+dbx*dbx+dby*dby)<<endl;
	  } else {
	    cout<<"WM(---)"<<endl;
	  }
	}
      }	
    }
    return thWB; 
  } 
  else {
    //if(inv_test_mode&&debugint){
    if(inv_test_mode){
      debug=true;
      cout<<"not wide"<<endl;
    }
    return invmapWittMao(p);
  }
};

///Inverse lense map in the wide-binary limit.
///
///When L>>1 (ie than Einstein angular radius) then we can pursue and interative solution in
///which the beta distance to (all-but) one of the lens objects can be considered far.  Which
///lens may be close is specified by (the sign of) iwhich.  Then the calculation proceeds by a Newton's method
///approach with the zeroth order solution being an exact single-lens solution with the distant
///lens object ignored.  Subsequent terms include an estimated contribution from the distant
///lens object. 
///
///Let \f$\vec p=p_{near}\f$ indicate the Einstein-scaled angular location of the possibly near lens object and
///define \f$\zeta=\beta_x-p_x + i(\beta_x-p_x)\f$
///We will solve for: \f$z=\theta_x -p_x + i(\theta_y-p_y)\f$ iteratevely. At each iteration we solve:
/// \f[
///    z\approx\zeta+\frac{\nu_{near}}{z^*}+\epsilon
/// \f]
/// where, \f$\epsilon\f$ represents the lens potential contribution from the distant lens object. 
/// For the zeroth iteration, we take \f$\epsilon=0\f$ and subsequently \f$\epsilon=\frac{\nu_{far}}{z+c}\f$
/// using the most recent estimate for \f$z\f$ and \f$c=p_{x,far}-p{x,near}=\pm L\f$.
/// At each iteration the solution is:
/// \f[
///    |z|\approx\frac{|\zeta+\epsilon|}2\left(\sqrt{1+\frac{4\nu_{near}}{|\zeta+\epsilon|^2}}\pm1\right)
/// \f]
/// with the complex argument given such that \f$ (\zeta+\epsilon)z^*\f$ is real.
/// As usual for a single lens there two image solutions, inside and outside the Einstein ring. 
/// These are selected above by the indicated choice of sign.  We iterate
/// separately a series of approximate solutions with either sign.
/// Where this approximation is appropriate, there will be one additional solution is near the other lens.
/// Note that it is straightforward to generalize to multiple lenses by changing the definition of \f$\epsilon\f$
/// to include contrinbutions from additional distant lenses.
///
vector<Point> GLensBinary::invmapWideBinary(const Point &p){
  const int maxIter=1000;
  //const int maxIter=20;
  //Lens equation:
  //double x=p.x,y=p.y,x1=x-L/2,x2=x+L/2,r1sq=x1*x1+y*y,r2sq=x2*x2+y*y;
  //double a,b,c1=(1-nu)/r1sq,c2=nu/r2sq;
  //return Point(x-x1*c1-x2*c2,y-y*(c1+c2));
  //pick dominant lens:
#ifdef USE_KIND_16
  typedef long double double_type;
#else
  typedef double double_type;
#endif
  double_type xL=p.x,yL=p.y,x1=xL-(double_type)L/2.0L,x2=xL+(double_type)L/2.0L,r1sq=x1*x1+yL*yL,r2sq=x2*x2+yL*yL;
  double_type nu_neg=(double_type)nu,nu_pos=1.0L-nu_neg,cpos=nu_pos*nu_pos/r1sq,cneg=nu_neg*nu_neg/r2sq;
  double_type nu_n,nu_f,c;
  vector<Point>result;
  if(cpos<cneg){//close (dominant) point is one x<0 half-plane
    //cout<<"neg-dominant"<<endl;
    c=L;
    nu_n=nu_neg;
    nu_f=nu_pos;
  } else {
    //cout<<"pos-dominant"<<endl;
    c=-(double_type)L;
    nu_n=nu_pos;
    nu_f=nu_neg;
  }
  complex<double_type> zp,zm,zf,zeta,ep,em,ef;
  zeta=complex<double_type>(p.x+c/2.0L,p.y);
  ep=em=ef=complex<double_type>(0,0);
  double err=1,relerr=1;
  int iter=0;
  //For initialization we consider two options.  Initialize with the single lens solution, the nominal approach described above,
  //or, when we already have a good solution, (A) initialize with the last-found solution.
  //  th[0]=Point(real(zp)-c/2.0L,imag(zp))   -->   zp = th[0].x + c/2 + I th[0].y
  //  th[1]=Point(real(zm)-c/2.0L,imag(zm))   -->   zp = th[1].x + c/2 + I th[1].y
  //  th[2](Point(real(zf)+c/2.0L,imag(zf))   -->   zp = th[2].x - c/2 + I th[2].y
  //  then,
  //  ep=ep(zp)
  //  em=em(zm)
  //  ef=ef(zf)
  //
  //  Or, alternatively, (B) ep=em=ef=0;
  // We run with the result with the smaller residual.  This also should avoid problems when there are changes in the root order,...
  //First the zero case, which is just one iteration of the main loop with ep=em=ef=0.
  //cout<<" Wide Binary x y = "<<xL<<" "<<yL<<"  c="<<c<<" nu="<<nu_n<<endl;
  if(save_thetas_wide){//This should be equivalent to the nominal result, but there are at least numerical differences.  Needs investigation...
    double p2c2=p.x*p.x+p.y*p.y+c*c;
    double pc=c*p.x;//note ym2=yp2=p2+pc,yf2=p2-pc;
    double_type root=sqrt(1.0L+4.0L*nu_n/(p2c2+pc));
    double_type rootf=sqrt(1.0L+4.0L*nu_f/(p2c2-pc));
    zp=zeta*(double_type)((1.0L+root)/2.0L);
    zm=zeta*(double_type)((1.0L-root)/2.0L);
    zf=(zeta-c)*(double_type)((1.0L-rootf)/2.0L);
    ep=nu_f/conj(zp-c);
    em=nu_f/conj(zm-c);
    ef=nu_n/conj(zf+c);
    //err=abs(ep)+abs(em)+abs(ef);
    err=abs(zp-zeta-nu_n/conj(zp)-ep);//comput residual of lens eq.
    err+=abs(zm-zeta-nu_n/conj(zm)-em);
    err+=abs(zf+c-zeta-nu_f/conj(zf)-ef);
    iter=1;
    //if(verbose) cout<<"save_thetas_wide p=*"<<p.x<<","<<p.y<<"):ep,em,ef"<<ep<<","<<em<<","<<ef<<endl;
  }
  if(save_thetas_wide and have_saved_soln and theta_save.size()==3){
    zp=complex<double_type>(theta_save[0].x+c/2.0L,theta_save[0].y);  //Need to change from saved_roots to theta_save...
    zm=complex<double_type>(theta_save[1].x+c/2.0L,theta_save[1].y);
    zf=complex<double_type>(theta_save[2].x-c/2.0L,theta_save[2].y);
    //we borrow the "old" variables for this.
    complex<double_type>epsave=nu_f/conj(zp-c);
    complex<double_type>emsave=nu_f/conj(zm-c);
    complex<double_type>efsave=nu_n/conj(zf+c);
    cout<<"          ep0="<<zp-zeta-nu_n/conj(zp)<<endl;
    cout<<"          em0="<<zm-zeta-nu_n/conj(zm)<<endl;
    cout<<"          ef0="<<zm+c-zeta-nu_f/conj(zf)<<endl;
    cout<<"  zp="<<zp<<"  ep="<<epsave<<endl;
    cout<<"  zm="<<zm<<"  em="<<emsave<<endl;
    cout<<"  zf="<<zf<<"  ef="<<efsave<<endl;
    double_type errsave=abs(zp-zeta-nu_n/conj(zp)-epsave);
    errsave+=abs(zm-zeta-nu_n/conj(zm)-emsave);
    errsave+=abs(zf+c-zeta-nu_f/conj(zf)-efsave);
    //And this is the test:
    cout<<" naive err="<<err<<"  save err="<<errsave<<endl;
    if(errsave<err){
      ep=epsave;
      em=emsave;
      ef=efsave;
      err=errsave;
    }
    cout<<"save_thetas_wide (chose):ep,em,ef"<<ep<<","<<em<<","<<ef<<endl;
    iter=1;
  }//now we have our preferred choice of ep/em/ef and we are ready for the main loop.

  bool fail=false;
  while(err>epsTOL&&relerr>epsTOL){
    //cout<<"WideBinary: iter="<<iter<<"\n          err="<<err<<"  relerr="<<relerr<<endl;
    iter++;
    if(iter>maxIter){
      if(debug)cout<<"invmapWideBinary maxIter reached: Failing."<<endl;
      return result;
    }
    //cout<<" err="<<err<<endl;
    complex<double_type>epold=ep,emold=em,efold=ef;
    complex<double_type>yp=zeta+ep,ym=zeta+em;
    complex<double_type> yf=zeta-c+ef;
    double_type yp2=real(yp*conj(yp));
    double_type ym2=real(ym*conj(ym));
    double_type yf2=real(yf*conj(yf));
    double_type zpfac=(sqrt(1.0L+4.0L*nu_n/yp2)+1.0L)/2.0L;
    double_type zmfac=-(sqrt(1.0L+4.0L*nu_n/ym2)-1.0L)/2.0L;
    double_type zffac=-(sqrt(1.0L+4.0L*nu_f/yf2)-1.0L)/2.0L;
    zp=yp*zpfac;
    zm=ym*zmfac;
    zf=yf*zffac;
    //if(debug){
    // cout<<"  zpfac="<<zpfac<<"  yp2="<<yp2<<" c="<<c<<endl;
    // cout<<"  zmfac="<<zmfac<<"  ym2="<<ym2<<" c="<<c<<endl;
    // cout<<"  zffac="<<zffac<<"  yf2="<<yf2<<" -c="<<-c<<endl;
    //}
    ep=nu_f/conj(zp-c);
    em=nu_f/conj(zm-c);
    ef=nu_n/conj(zf+c);
    err=abs(ep-epold)+abs(em-emold)+abs(ef-efold);
    double_type zmean=abs(zp)+abs(zm)+abs(zf);
    relerr=abs((ep-epold)*zp*zp/zmean/zmean)+abs((em-emold)*zm*zm/zmean/zmean)+abs((ef-efold)*zf*zf/zmean/zmean);
    if(debug){
      cout<<"  zp="<<zp<<"  ep="<<ep<<" D="<<ep-epold<<endl;
      cout<<"  zm="<<zm<<"  em="<<em<<" D="<<em-emold<<endl;
      cout<<"  zf="<<zf<<"  ef="<<ef<<" D="<<ef-efold<<endl;
      cout<<" err,relerr="<<err<<","<<relerr<<" vs "<<epsTOL<<endl;
      cout<<"  mp="<<mag(Point(real(zp)-c/2.0L,imag(zp)))
	  <<"  mm="<<mag(Point(real(zm)-c/2.0L,imag(zm)))
	  <<"  mf="<<mag(Point(real(zf)-c/2.0L,imag(zf)))<<endl;
    }
  }
  //Should we order the roots for consistency with WittMao results??? 
  result.push_back(Point(real(zp)-c/2.0L,imag(zp)));
  result.push_back(Point(real(zm)-c/2.0L,imag(zm)));
  result.push_back(Point(real(zf)+c/2.0L,imag(zf)));
  //perhaps save the result.
  if(save_thetas_wide ){
    theta_save=result;
    have_saved_soln=true;
  } else {
    have_saved_soln=false;
  }
  return result;
};

vector<Point> GLensBinary::invmapWittMao(const Point &p,bool no_check){
  bool test_images=true;
  //if no_check==true then we return all polynomial roots without checking that they are consistent with the forward map.
  //Here quantities are like WittMao95, but with:
  // z         ->  z*z1
  // 2m        -> M*z1^2  ; but in our normalization 2 m = nu+(1-nu) = 1
  // 2 Delta m -> DM*z1^2 ; in our normalization 2 Delta m = nu - (1-nu) = 2*nu-1
  // zeta      -> B*z1
  const complex<double> cplx_i=complex<double>(0,1),cplx_1=1;
  double z1=L/2;
  double z12=z1*z1;
  double M=1/z12;
  double DM=M*(2*nu-1);
  complex<double> B=(p.x+cplx_i*p.y)/z1,Bb=conj(B);
  complex<double> u,v,w,wb;
  complex<double> c[7],roots[6];
  int nroots=5;
  bool z1scaled=false;
  //if(z1>LEADTOL){    //scaled with z->z*z1,ck->ck*z1^-k
  if(z1>1){
    z1scaled=true;//scaled with z->z*z1,ck->ck*z1^-k
    u=DM*B+M;
    w=M*Bb+DM;
    wb=conj(w);
    v=wb+Bb;//M*B+DM+Bb
    c[5]=1.0-Bb*Bb;
    c[4]=-B*c[5]-w;
    c[3]=2.0*v*Bb-2.0;
    c[2]=M*wb-2.0*(Bb*u+c[4]);
    c[1]=1.0-v*v-M*M*conj(c[5]);
    c[0]=(DM+2.0*Bb)*u+c[4];
  } 
  else if(true) {    //scaled with z->z*z1,ck->ck*z1^-k
    B=(p.x+cplx_i*p.y);Bb=conj(B);   //=B1*z1     
    M=1;
    DM=2*nu-1;            //=(DM*z1^2 above)
    u=DM*B+z1;           //(=u*z1^3 above)
    w=Bb+DM*z1;        //(=w*z1^3 above)
    wb=conj(w);
    v=wb+Bb*z12;         //(=v*z1^3 above)
    c[5]=z12-Bb*Bb;                         //c5*z1^2
    c[4]=-B*c[5]-w;                         //c4*z1^3
    c[3]=2.0*v*Bb-2.0*z12*z12;              //c3*z1^4
    c[2]=wb-2.0*(Bb*u*z1+c[4]*z12);         //c2*z1^5
    c[1]=z12*z12*z12-v*v-conj(c[5]);        //c1*z1^6
    c[0]=((DM+2.0*Bb*z1)*u+c[4]*z12)*z12;//c0*z1^7
  } 
  else { //handle degenerate (single point lens) case without z1 scalingelse { //handle degenerate (single point lens) case without z1 scaling
    B=(p.x+cplx_i*p.y);Bb=conj(B);
    M=1;
    DM=2*nu-1;
    c[5]=0;c[4]=0;
    c[3]=-Bb*Bb;
    c[2]=-B*c[3]-Bb;
    c[1]=2.0*B*Bb;
    c[0]=B;
    nroots=3;
  }
  //Test for effectively lower order; 
  //If the leading order polynomial coefficients are effectively zero, then the solve will fail.
  double c_lower=0;
  for(int k=0;k<nroots;k++)c_lower+=abs(real(c[k]))+abs(imag(c[k]));
  for(int i=nroots;i>1;i--){
    double val=(abs(real(c[i]))+abs(imag(c[i])))/c_lower;
    if(debug)cout<<"val["<<i<<"]= "<<val<<" ?> "<<LEADTOL<<endl;
    if(val>LEADTOL)break;
    //break; //debug  Seems to work better without throwing out terms here  WHY?
    //The original reason for this was failure of the polynomial solution when the
    //leading order coeficient vanishes.
    //I don't think I know how it behaves when the coeficient is just very small, though it makes sense that
    //there must at least be an underflow issue at some point.
    //I think, at the time, I wasn't testing regions of parameter space with such extreme values..
    nroots--;
    c_lower-=abs(real(c[nroots]))+abs(imag(c[nroots]));
  }
  //Debug
  if(debug){
    cout<<"nu="<<nu<<" z1="<<z1<<" M="<<M<<" DM="<<DM<<" B= "<<B<<" u="<<u<<" w="<<w<<" v="<<v<<endl;
    cout<<"nroots="<<nroots<<endl;
    for( int i=0;i<=5;i++)
      cout<<"c["<<i<<"]="<<c[i]<<(i>nroots?"**":"")<<endl;
  }

  bool roots_changed =true;//(this isn't used but compiler may be concerned about an uninitialized value) 
  bool polish_only = false;
  if(nroots==5 and save_thetas_poly and have_saved_soln and theta_save.size()==5){
    //for(int i=0;i<5;i++)roots[i]=saved_roots[i];
    for(int i=0;i<5;i++)roots[i]=complex<double>(theta_save[i].x,theta_save[i].y);
    polish_only=true;
    if(verbose){
      cout<<"using saved roots: p=("<<p.x<<","<<p.y<<")"<<endl;
      for(int i=0;i<5;i++)cout<<"  "<<i<<": "<<roots[i]<<endl;
    }
  } else if(save_thetas_poly and verbose) cout<<"   no saved roots: p=("<<p.x<<","<<p.y<<")"<<endl;

  if(nroots==5)cmplx_roots_5(roots, roots_changed, c, polish_only);
  else cmplx_roots_gen(roots, c, nroots,true,false);
  if(save_thetas_poly and nroots==5){
    theta_save.resize(5);
    for(int i=0;i<5;i++)theta_save[i]=Point(real(roots[i]),imag(roots[i]));
    have_saved_soln=true;
  } else {
    have_saved_soln=false;
  }
  
  if(nroots==4&&fix_nr4_roots){//try to make linear correction to root values
    for( int i=0;i<nroots;i++){
      complex<__float128> dP4=0.0;
      complex<__float128> z=roots[i];
      for(int k=1;k<=4;k++){
	dP4+=c[k]*(double)k;
	dP4/=z;
      };
      //Result is dP4 = dP/dz(z_k) / z_k^4
      //dP4=((((c1/z+2*c2)/z+3*c3)/z+4*c4)/z
      //   =(c1+2*c2*z+3*c3*z^2+4*z4+z^4)/z^4
      complex<double> double_dP4;double_dP4=dP4;
      complex<double> delta=-c[5]/double_dP4;
      if(debug)cout<<"dP4="<<double_dP4<<" delta="<<delta<<endl;
      ///We can derive this following expression by setting P5(x+eps)-P4(x)=0, then solve liniarly for eps
      ///For abs(delta)<<0.2 and abs(delta)>>0.2 the result approximates that in the comment above, but all values other than
      ///delta=-0.2 are allowed.  In this case, the 0.2 arises because n=5, not arbitrarily.
      roots[i]*=(1.0+cplx_1/(5.0+cplx_1/delta));
    }
  }
  //For debugging, test the results;
  if(debug){
    cout<<"complex poly test:"<<endl;
    for( int i=0;i<nroots;i++){
      complex<__float128> res=0.0;
      complex<__float128> z=roots[i];
      for(int k=5;k>=0;k--){
	res*=z;
	res+=c[k];
      }
      complex<long double> lres,lz;
      lres=res;lz=z;
      cout<<i<<": "<<lz<<"->"<<lres<<endl;
    }
  }
  
  vector<Point> result;
  const double TOL=LEADTOL*LEADTOL;
  double maxerr=0;
  int imaxerr=-1;
  double minerrfail=INFINITY;
  int iminerrfail=-1;

  for(int i=0;i<nroots;i++){
    Point newp=Point(0,0);
    if(z1scaled)
      newp=Point(real(roots[i])*z1,imag(roots[i])*z1);
    else 
      newp=Point(real(roots[i]),imag(roots[i]));
    if(no_check) result.push_back(newp);
    else{
      //Test solutions: Could also try keeping those with err/(1+B*Bb)<TOL say
      Point btheta=map(newp);
      double dx=btheta.x-p.x,dy=btheta.y-p.y;
      double err=dx*dx+dy*dy;
      //double dbx=newp.x-p.x,dby=newp.y-p.y;
      double x=newp.x,y=newp.y,x1=x-L/2,x2=x+L/2,r1sq=x1*x1+y*y,r2sq=x2*x2+y*y;
      double c1=(1-nu)/r1sq,c2=nu/r2sq;
      err/=(1+c1+c2);
      if(debug ){
	cout<<"testing: "<<roots[i]<<" -> |("<<btheta.x<<"-"<<p.x<<","<<btheta.y<<"-"<<p.y<<")|="<<sqrt(dx*dx+dy*dy)<<" "<<(dx*dx+dy*dy<TOL*(1+c1+c2)?" < ":"!< ")<<LEADTOL*sqrt(1+c1+c2)<<"="<<LEADTOL<<"*√(1+"<<c1<<"+"<<c2<<")"<<endl;
      }
      if(err<TOL){     //RHS is squared estimate in propagating error of LEADTOL in root through map()
	result.push_back(newp);
	if(err>maxerr){
	  maxerr=err;
	  imaxerr=result.size()-1;
	}
      } else {
	if(err<minerrfail){
	  minerrfail=err;
	  iminerrfail=i;
	}
      }
      //For now we just adopt the ordering from the SG code.  Might change to something else if needed...
    };
  }
  if(test_images){
    ///After first pass application of the lens map test to each images we apply checks on the
    ///overall set of images.  There should be at least NimageMin and no fewer than NimageMax
    ///images and the count of any images in excess of NimageMin should be even.
    ///If there are too many images, then we discard the one which most-marginally passed the map
    ///test.  If there is one too few, then we add back a candidate image which moderately failed the map test (within a factor of 100), if present.  We could potentially do better by checking
    ///the parities of the images to be added or deleted.
    ///Note, that later, for finite source cases, we may add back more images deemed to be
    ///missing extremely near the lens points.  That logic could be applied here as well.
    int ni=result.size();
    bool pass1=( ni >= NimageMin );
    bool pass2=(ni <= NimageMax);
    bool pass3=( (ni-NimageMin)%2==0 );
    if(not (pass1 and pass2 and pass3)){
      if(pass1 and not pass3){//This is something we can try to fix:
	//cout<<"erasing a candidate image "<<imaxerr<<endl;
	result.erase(result.begin()+imaxerr);
      } else if(ni+1==NimageMin and minerrfail<TOL*100){//Looks like we missed an image near one of the lenses, but can see it with related TOL
	Point newp=Point(0,0);
	int i=iminerrfail;
	if(z1scaled)
	  newp=Point(real(roots[i])*z1,imag(roots[i])*z1);
	else 
	  newp=Point(real(roots[i]),imag(roots[i]));
	double mg=mag(newp);
	if(mg<0)result.push_back(newp);
      }else {
	cout<<"something wrong with the image set "<<pass1<<" "<<pass2<<" "<<pass3<<", but we forge on"<<endl;
      }
    }
  }

  if(debug)cout<<"Found "<<result.size()<<" images."<<endl;

  return result;
};

int GLensBinary::poly_root_integration_func_vec (double t, const double theta[], double thetadot[], void *instance){
  //This is similar to the image point integration func, but we integrate the roots of the poly to be directly
  //We use notation beta.x + i beta.y -> w and theta.x + i theta.y -> z
  //
  //Another change is that we have here written this function explicitly for sepcific derived GLens class
  //
  GLensBinary *thisobj = static_cast<GLensBinary *>(instance);
  const Trajectory *traj=thisobj->trajectory;
  Point beta0=thisobj->get_obs_pos(*traj,t);
  Point betadot=thisobj->get_obs_vel(*traj,t);
  complex<double> wdot(betadot.x,betadot.y);
  bool fail=false;
  for(int k=0;k<2*thisobj->NimageMax;k++)thetadot[k]=0;
  //cout<<"beta0=("<<beta0.x<<","<<beta0.y<<")"<<endl;
  for(int image =0; image<NimageMax;image++){
    //Point p(theta[image*2+0],theta[image*2+1]);
    complex<double> z(theta[image*2+0],theta[image*2+1]);
    //double j00i,j10i,j01i,j11i,invJ = thisobj->invjac(p,j00i,j01i,j10i,j11i);
    /*
    //This block was provides contraint-driving in the original formulation
    //It seems possible, but perhaps complicated to do this in the poly_root case
    //The issue is that there is no explict inverse transform w(z) with the polynomial
    //because that polynomial includes up to cubic-order combinations of w and w*.
    //Using the original map equation w(z), as before, is not appropriate here since
    //The polynomial has solutions which are inconsistent with that map function.
    //A plausible way to implement this is to assume that the error in w is small
    //then solve the w(z),w*(z) system in the small dw limit.  The resulting expression
    //is a little bit complicated, so we don't yet implement that.
    Point beta=thisobj->map(p);
    double dx=beta.x-beta0.x,dy=beta.y-beta0.y;
    double adjbetadot[2];
    //double rscale=thisobj->estimate_scale(Point(theta[0],theta[1]));
    //if our image is near one of the lenses, then we need to tread carefully
    double rk=thisobj->kappa;//*rscale;
    adjbetadot[0] = betadot.x;
    adjbetadot[1] = betadot.y;
    if(isfinite(dx))adjbetadot[0] += -rk*dx;
    if(isfinite(dy))adjbetadot[1] += -rk*dy;
    */
    //Next we compute derivative, which previously looked like:
    //
    //if(isfinite(invJ)){
    //thetadot[2*image]   = (j00i*adjbetadot[0]+j01i*adjbetadot[1]);
    //thetadot[2*image+1] = (j10i*adjbetadot[0]+j11i*adjbetadot[1]);
    //}
    //
    //In the complex notation, the effective invjacobian is simpler
    // zdot = ( wdot + Ecc*wdotcc) * mu
    // where E = nu1/z1^2 + nu2/z2^2, cc indicates conjugate and
    // the magnification mu = 1 / (1 - E*Ecc) 
    complex<double> c1((thisobj->L)/2,0), z1=z-c1, z2=z+c1;
    complex<double> E = (1-nu)/z1/z1 + nu/z2/z2;
    complex<double> zdotinvmu =  wdot +conj(E*wdot) ;
    double invmu = 1-norm(E);
    if(invmu==0){
      //here we panic and rather than return NAN, we evolve as if the lensing effect is unit scale
      //this is wrong, but it will move use past the point where
      invmu=1;
    }
    thetadot[2*image]   = zdotinvmu.real()/invmu;
    thetadot[2*image+1] = zdotinvmu.imag()/invmu;
  }

  return !fail?GSL_SUCCESS:3210124;
}

double GLensBinary::mag(const Point &p){
  double x=p.x,y=p.y,x1=x-L/2,x2=x+L/2,r1sq=x1*x1+y*y,r2sq=x2*x2+y*y;
  //cout<<"x,y,r1sq,r2sq:"<<x<<" "<<y<<" "<<r1sq<<" "<<r2sq<<endl;
  double m1=(1-nu), m2=nu, dEr4=m1*r2sq-m2*r1sq;
  double cosr2=x1*x2+y*y,cos2r4=cosr2*cosr2;
  double r4=r1sq*r2sq,r8=r4*r4;
  double invmuR8=r8-dEr4*dEr4-4*m1*m2*cos2r4;
  double mg=0;
  if(r8>0)mg=r8/invmuR8;
  //cout<<"mg="<<mg<<endl;
  return mg;
};

///returns J=det(d(map(p))/dp)^-1, sets, j_ik = d(map(pi))/dpk
double GLensBinary::jac(const Point &p,double &j00,double &j01,double &j10,double &j11){
  double x=p.x,y=p.y,x1=x-L/2,x2=x+L/2,y2=y*y,r1sq=x1*x1+y2,r2sq=x2*x2+y2;
  double E1=(1-nu)/r1sq, E2=nu/r2sq,dE=E1-E2,E=E1+E2;
  double cos2=x1*x2+y2;cos2*=cos2/r1sq/r2sq;
  double invmu=1-dE*dE-4*E1*E2*cos2;
  double mu=1/invmu;
  double E1w=2*E1/r1sq,E2w=2*E2/r2sq;
  double f=(E1w+E2w)*y2 - E;
  j10 = j01 = (E1w*x1+E2w*x2)*y;
  j00       = 1-f;
  j11       = 1+f;
  return mu;
};

///The inverse jacobian can be finite when the jacobian is not so we compute it directly
double GLensBinary::invjac(const Point &p,double &ij00,double &ij01,double &ij10,double &ij11){
  double x=p.x,y=p.y,x1=x-L/2,x2=x+L/2,y2=y*y,r1sq=x1*x1+y2,r2sq=x2*x2+y2;
  double E1r4=(1-nu)*r2sq, E2r4=nu*r1sq,dEr4=E1r4-E2r4,Er4=E1r4+E2r4;
  double cosr2=x1*x2+y2,cos2r4=cosr2*cosr2;
  double r4=r1sq*r2sq,r8=r4*r4;
  double invmuR8=r8-dEr4*dEr4-4*nu*(1-nu)*cos2r4;
  double mu=r8/invmuR8;
  double E1wr8=2*E1r4*r2sq,E2wr8=2*E2r4*r1sq;
  double fr8=(E1wr8+E2wr8)*y2 - Er4*r4;
  ij10 = ij01 = -(E1wr8*x1+E2wr8*x2)*y/invmuR8;
  ij00       = (r8+fr8)/invmuR8;
  ij11       = (r8-fr8)/invmuR8;
  return mu;
};

///Compute the complex lens shear, and some number of its derivatives
/// gamma = \sum_i^N nu_i / (zc*zc)   =  dbetac/dz
vector<complex<double> > GLensBinary::compute_shear(const Point &p, int nder)const{
  double x=p.x,y=p.y,x1=x-L/2,x2=x+L/2;
  complex<double> z1(x1,y), z2(x2,y);
  vector<complex<double> >gammas;   
  complex<double> dNgamma1=(1-nu)/z1/z1, dNgamma2=nu/z2/z2;
  gammas.push_back(dNgamma1+dNgamma2);
  for(int n=0;n<nder;n++){
    dNgamma1 *=-(n+2.0)/z1;
    dNgamma2 *=-(n+2.0)/z2;
    gammas.push_back(dNgamma1+dNgamma2);
  }
  return gammas;
}
 
///Compute the Laplacian of the local image magnification explicitly
///
///See explanation with GLens::Laplacian_mu
/*
double GLensBinary::Laplacian_mu(const Point &p)const{
  //The first part is specialized to the binary lens, the rest is general...
  double x=p.x,y=p.y,x1=x-L/2,x2=x+L/2,y2=y*y,r1sq=x1*x1+y2,r2sq=x2*x2+y2;
  complex<double> z1(x1,y1), z2(x2,y2);
  complex<double>   gamma1=(1-nu)/z1/z1,    gamma2=nu/z2/z2,       gamma=gamma1+gamma2, gammac=conj(gamma);
  complex<double>  dgamma1=-2*gamma1/z1,   dgamma2=-2*gamma2/z2,  dgamma=dgamma1+dgamma2;
  complex<double> d2gamma1=-3*dgamma1/z1, d2gamma2=-3*gamma2/z2, d2gamma=d2gamma1+d2gamma2;

  //The rest is general for any lens system
  double invmu = 1-norm(gamma);
  if(invmu==0)return 0;//Fail gracefully;  There is no appropriate result near caustics
  double mu=1/invmu, mu2=mu*mu;
  double term1 = real(-2*gammac*gammac*( 3*mu2*gammac*dgamma*dgamma + mu*d2gamma ));
  double term2 = real(norm(dgamma)*( 6*mu*(mu-1) + 1 ));
  return 4*mu*mu2*(term1+term2);
};
*/

