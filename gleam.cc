//Gravitational lens event analysis machine
//Written by John G Baker NASA-GSFC (2014)

//#include "mlfit.hh"
#include <valarray>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include "omp.h"
#include "options.hh"
//#include <mcheck.h>
#include "bayesian.hh"
#include "proposal_distribution.hh"
#include "ptmcmc.hh"
#include "mldata.hh"
#include "mlsignal.hh"
#include "mllike.hh"

using namespace std;

typedef initializer_list<double> dlist;
typedef initializer_list<int> ilist;

bool debug = false;
bool debugint = false;
bool debug_signal = false;

shared_ptr<Random> globalRNG;//used for some debugging... 

//Global Control parameters (see more in main())
const int Npar=9;
const double MaxAdditiveNoiseMag=22;
int output_precision;
double mm_lens_rWB;

//Analysis functions defined below.
void dump_view(const string &outname, bayes_data &data, ML_photometry_signal &signal, bayes_likelihood &like,state &s,double tstart,double tend,int nsamples);
void dump_mag_map(const string &outname,bayes_data &data,ML_photometry_signal &signal, state &s,double tstart,double tend,int nsamples=301,int cent=-2,bool output_nlens=false);
void dump_lightcurve(const string &outname,bayes_likelihood&like,state &s,double tstart,double tend,int nsamples=301);


//***************************************************************************************8
//main test program
int main(int argc, char*argv[]){
  //string datafile;
  const int NparRead=Npar; 

  //Create the sampler
  ptmcmc_sampler mcmc;
  bayes_sampler *s0=&mcmc;
  //Create the model components
  GLensBinary binarylens;
  Trajectory linear_trajectory(Point(0,0), Point(1,0));
  Trajectory *traj=&linear_trajectory;
  GLens *lens=&binarylens;

  ///We are at the point now that we need to generalize the data types on the command line
  ///So far there are 2 tested types ML_OGLE_data and ML_generic_data while ML_photometry_mock_data is
  ///under development.  My plan is to define a parent-class pointer ML_photometry *data, and to 
  ///then instantiate it based on what was passed in.  This means that we will have to parse options 
  ///twice.  First just to see what type of data we have, with the first call to parse. We can ignore
  ///any failure notice at that point.  Next we add the appropriate options the Options object before
  ///parsing again.  This will mean that we need a common "setup" interface, with no data-type specific
  ///arguments.  Any such arguments will need to come from the command line.  I was thinking the
  ///options may be such as -OGLEdata=<filename> or -mock_data (with mock_tstart=..., etc for other params).
  ///This same command-line interface can later be used for joint anaysis of multiple data sets.  All we
  ///will need for that may be a compounding type bayes_data class.  [Thinking ahead, we may later have
  ///astrometric data as well. These might be similarly added, but they will involve a different type of
  ///signal model.  We may need to be able to use vector-valued bayes_data::values with a vector of type
  ///channel-labels that can be checked.  Maybe the channels are just strings, but if there is important 
  ///meta-information regarding the channels, then these could become instances of a channel class.]
  Options opt;
  ML_photometry_data::addStaticOptions(opt);
  //For the first-pass option check, we first make a copy of argc/argv
  int ac=argc;char* av[argc];
  char * filename=NULL;
  for(int i=0;i<argc;i++){av[i]=argv[i];};
  opt.parse(ac,av,false);
  //now select the data obj
  ML_photometry_data *data;
  if(opt.set("OGLE_data"))
    data=new ML_OGLEdata();
  else if(opt.set("gen_data"))
    data=new ML_generic_data();
  else if(opt.set("mock_data"))
    data=new ML_mock_data();
  else {
    //for backward compatibility [deprecated] default is to assume OGLE data and try to read the data from a file named in the (extra) first argument
    if(argc>=1){
      cout<<"Setting filename from first argument for backward compatibility [deprecated]."<<endl;
      filename=av[1];
      data=new ML_OGLEdata;
    } else { //go on assuming mock data; Except for backward compatibility, this would be the default.
      cout<<"No data file indicated!"<<endl;
      data=new ML_mock_data();
    }
  }
      
  //Eventually want to handle signal polymorphism similarly
  ML_photometry_signal signal(traj, lens);

  s0->addOptions(opt,"");
  lens->addOptions(opt,"");
  data->addOptions(opt,"");
  signal.addOptions(opt,"");

  opt.add(Option("nchains","Number of consequtive chain runs. Default 1","1"));
  opt.add(Option("seed","Pseudo random number grenerator seed in [0,1). (Default=-1, use clock to seed.)","-1"));

  //prior options:
  opt.add(Option("log_tE","Use log10 based variable (and Gaussian prior with 1-sigma range [0:log10(tE_max)] ) for tE parameter rather than direct tE value."));
  opt.add(Option("tE_max","Uniform prior max in tE. Default=100.0/","100.0"));

  //likelihood or data option?
  opt.add(Option("additive_noise","Interpret Fn as magnitude of additive noise. Fn_max is magnitude of maximum noise level (i.e. minimum noise magnitude)"));
  opt.add(Option("Fn_max","Uniform prior max (min for additive) in Fn. Default=1.0 (18.0 additive)/","1"));
  opt.add(Option("view","Don't run any chains, instead take a set of parameters and produce a set of reports about that lens model."));

  //other options
  opt.add(Option("magmap","Don't run any chains, instead just make a magnitude map.  In this case the parameters should be just q,L,width."));
  opt.add(Option("mm_center","On which lens to center magmap. (-1,0,1), with default zero for CoM.","0"));
  opt.add(Option("mm_d0x","Explicit x coord offset for magmap center, with default zero.","0"));
  opt.add(Option("mm_d0y","Explicit y coord offset for magmap center, with default zero.","0"));
  opt.add(Option("mm_samples","Number of samples in magmap (default 300)","300"));
  opt.add(Option("mm_nimage","Include number of images magmap"));
  opt.add(Option("precision","Set output precision digits. (Default 13).","13"));
  opt.add(Option("mm_lens_rWB","In magmap mode, radial scale outside which we use the WideBinary lens inversion (neg for none). Use for quality control. (Default 5).","5"));

  int Nlead_args=1;
  if(filename)Nlead_args++;

  bool parseBAD=opt.parse(argc,argv);
  if(parseBAD||(argc != Nlead_args+1 && argc!=Nlead_args+1+Npar && argc!=5)) {
    cout << "You gave " << argc-1 << " arguments. Expecting "<<Nlead_args<<" or "<<Nlead_args+Npar<<" or 4."<< endl;
    cout << "Usage:\n gleam [-options=vals] data_file_name output_name [ I0 Fs Fn logq logL r0 phi tE tpass ]" << endl;
    cout << "Or:\n gleam -magmap [-options=vals] output_name logq logL width" << endl;
    cout <<opt.print_usage()<<endl;
    return 1;
  }
    
  cout<<"flags=\n"<<opt.report()<<endl;

  //Post parse setup
  lens->setup();  
  double nburn_frac,Tmax,Fn_max,tE_max,tcut,seed;
  int Nchain;
  int mm_center,mm_samples,save_every;
  double mm_d0x,mm_d0y;
  bool do_magmap,view;
  int Nsigma=1;
  int Nbest=10;
  view=opt.set("view");
  istringstream(opt.value("nchains"))>>Nchain;
  istringstream(opt.value("seed"))>>seed;
  //if seed<0 set seed from clock
  if(seed<0)seed=time(NULL);
  istringstream(opt.value("tcut"))>>tcut;
  do_magmap=opt.set("magmap");
  istringstream(opt.value("mm_center"))>>mm_center;
  istringstream(opt.value("mm_d0x"))>>mm_d0x;
  istringstream(opt.value("mm_d0y"))>>mm_d0y;
  istringstream(opt.value("mm_samples"))>>mm_samples;
  istringstream(opt.value("precision"))>>output_precision;
  istringstream(opt.value("mm_lens_rWB"))>>mm_lens_rWB;
  //Prior params (should move out to objects which manage specific parameters.)
  istringstream(opt.value("Fn_max"))>>Fn_max;
  istringstream(opt.value("tE_max"))>>tE_max;

  //read args
  string outname;
  ostringstream ss("");
  valarray<double>params(Npar);
  bool have_pars0=false;
  if(do_magmap&&argc==5){//different parameters in this case
    //datafile="";
    outname=argv[1];
    if(Npar<3){
      cout<<"We expected Npar>=3"<<endl;
      exit(1);
    }
    for(int i=0;i<3;i++)stringstream(argv[i+2])>>params[i];  
    have_pars0=true;
  } else if(argc>=Nlead_args+1){
    outname=argv[Nlead_args];
    if(argc>Nlead_args+1){
      for(int i=0;i<NparRead;i++)stringstream(argv[i+1+Nlead_args])>>params[i];
      have_pars0=true;
    }
  }

  //report
  cout.precision(output_precision);
  cout<<"\noutname = '"<<outname<<"'"<<endl;
  cout<<"seed="<<seed<<endl; 
  cout<<"Running on "<<omp_get_max_threads()<<" thread"<<(omp_get_max_threads()>1?"s":"")<<"."<<endl;

  //Should probably move this to ptmcmc/bayesian
  ProbabilityDist::setSeed(seed);
  globalRNG.reset(ProbabilityDist::getPRNG());//just for safety to keep us from deleting main RNG in debugging.
  
  //Set up ML objects
  signal.setup();
  //special handling for backward compatibility [deprecated]
  if(filename&&!do_magmap)dynamic_cast< ML_OGLEdata* >(data)->setup(filename);
  else data->setup();

  //Create the data object (This block is setting up the model.  If/when we separate model from data, then data obj may be defined below)
  double q0;
  istringstream(opt.value("q0"))>>q0;
  bool use_remapped_r0,use_remapped_q,use_log_tE,use_additive_noise=false;
  use_remapped_r0=opt.set("remap_r0");
  use_remapped_q=opt.set("remap_q");
  use_additive_noise=opt.set("additive_noise");
  use_log_tE=opt.set("log_tE");

  double r0s=6.0;
  if(use_remapped_r0){
    r0s=1.0;
  }

  //Set up the parameter space
  stateSpace space(Npar);
  {
    string names[]={"I0","Fs","Fn","logq","logL","r0","phi","tE","tpass"};
    if(use_additive_noise)names[2]="Mn";
    if(use_remapped_r0)names[5]="s(r0)";
    if(use_remapped_q)names[3]="s(1+q)";    
    if(use_log_tE)names[7]="log(tE)";
    space.set_names(names);  
    space.set_bound(6,boundary(boundary::wrap,boundary::wrap,0,2*M_PI));//set 2-pi-wrapped space for phi.
  }
  cout<<"&space="<<&space<<endl; 
  cout<<"Parameter space:\n"<<space.show()<<endl;

  //Handle separate magmap (only) option.
  if(do_magmap){//translate parameters
    double q,L,width;
    q=pow(10.0,params[0]);
    L=pow(10.0,params[1]);
    width=params[2];
    stateSpace lensSpace=lens->getObjectStateSpace();
    lens->defWorkingStateSpace(lensSpace);
    state lens_state(&lensSpace,valarray<double>({q,L}));
    lens->setState(lens_state);
    double x0=lens->getCenter(mm_center).x;
    cout<<"cent="<<mm_center<<" x0-xcm="<<x0+mm_d0x<<" xcm="<<lens->getCenter().x<<endl;
    Point pstart(x0+mm_d0x-width/2,mm_d0y-width/2);
    Point pend(x0+mm_d0x+width/2,mm_d0y+width/2);
    {
      ss.str("");ss<<outname<<"_mmap.dat";
      ofstream out(ss.str());
      out.precision(output_precision);
      lens->writeMagMap(out, pstart, pend, mm_samples);
    }
    if(opt.set("mm_nimage")){
      ss.str("");ss<<outname<<"_nmap.dat";    
      ofstream out(ss.str()); 
      out.precision(output_precision);
      lens->writeMagMap(out, pstart, pend, mm_samples,true);
    }
    exit(0);
  }    

  //Prune data
  cout<<"Ndata="<<data->size()<<endl;
  double t0,twidth;
  double tstart,tend;
  t0=data->getFocusLabel();
  data->getDomainLimits(tstart,tend);
  twidth=300;//Changed for study, may return to smaller range...;
  //twidth=10;
  //cout<<"ts="<<tstart<<" < t0="<<t0<<" < te="<<tend<<endl;

  //Initial parameter state
  state instate(&space,params);
  if(have_pars0){
    cout<<"Input parameters:"<<instate.show()<<endl;
  }


  //Set the prior:
  //Eventually this should move to the relevant constitutent code elements where the params are given meaning.
  const int uni=mixed_dist_product::uniform, gauss=mixed_dist_product::gaussian, pol=mixed_dist_product::polar; 
  //                                     I0      Fs     Fn          logq*   logL       r0     phi      tE     tpass  
  valarray<double>    centers((dlist){ 18.0,    0.5,   0.5*Fn_max,   0.0,   0.0,  r0s/2.0,   M_PI, tE_max/2,  t0     });
  //valarray<double> halfwidths((dlist){  5.0,    0.5,   0.5*Fn_max,   4.0,   2.0,  r0s/2.0,   M_PI, tE_max/2,  twidth });
  valarray<double> halfwidths((dlist){  5.0,    0.5,   0.5*Fn_max,   4.0,   1.0,  r0s/2.0,   M_PI, tE_max/2,  twidth });
  valarray<int>         types((ilist){gauss,    uni,   uni,          uni, gauss,      uni,    uni,      uni,  gauss  });
  if(use_remapped_q){//* 
    double qq=2.0/(q0+1.0);
    double ds=0.5/(1.0+qq*qq); //ds=(1-s(q=1))/2
    centers[3]=1.0-ds;
    halfwidths[3]=ds;          //ie range=[s(q=1),s(q=inf)=1.0]
    types[3]=uni;
  }
  if(use_log_tE){
    centers[7]=log10(tE_max)/2;
    halfwidths[7]=log10(tE_max)/2;
    types[7]=gauss;
  }
  if(use_additive_noise){
    if(Fn_max<=1)Fn_max=18.0;
    double hw=(MaxAdditiveNoiseMag-Fn_max)/2.0;
    centers[2]=MaxAdditiveNoiseMag-hw;
    halfwidths[2]=hw;
  }
  //**test HACK **// //Here we set the q/L params to near particular values for testing with Jeremy's data
  //centers[3]=log10(0.9);
  //halfwidths[3]=0.05;
  //types[3]=uni;
  //centers[4]=log10(0.6);
  //halfwidths[4]=0.05;
  //types[4]=uni;  
  //**end HACK **//

  mixed_dist_product prior(&space,types,centers,halfwidths);
  cout<<"Prior is:\n"<<prior.show()<<endl;

  //Set the likelihood
  //Should some of this move up before addOptions?
  bayes_likelihood *llike=nullptr;
  bayes_likelihood *like=nullptr;
  ML_photometry_likelihood mpl(&space, data, &signal, &prior);
  mpl.addOptions(opt,"");
  mpl.setup();
  like = &mpl;
  like->defWorkingStateSpace(space);

  ///At this point we are ready for analysis in the case that we are asked to view a model
  ///Note that we still have needed the data file to create the OGLEdata object, and concretely
  ///to set the domain.  This could be changed...
  if(view){
    if(have_pars0){
      int nsamples;
      double tfinestart,tfineend;
      like->getFineGrid(nsamples,tfinestart,tfineend);
      cout<<"Producing report on the model with specified parameters."<<endl;
      dump_view(outname, *data, signal, mpl, instate, tfinestart, tfineend, mm_samples);
    } else {
      cout<<" The -view option requires that parameters are provided."<<endl;
    }
  }

  //Just a little reporting
  if(have_pars0){
    double ll=like->evaluate_log(instate);
    double lp=prior.evaluate_log(instate);
    cout<<"log-Likelihood at input parameters = "<<ll<<endl;
    cout<<"log-posterior at input parameters = "<<ll+lp<<endl;
  }
  if(view)exit(0);
	  
  //assuming mcmc:
  //Set the proposal distribution (might want to formalize some of this, now basically copied in similar form in the for the third time in a new program...)
  int Ninit;
  proposal_distribution *prop=ptmcmc_sampler::new_proposal_distribution(Npar,Ninit,opt,&prior,&halfwidths);
  cout<<"Proposal distribution (at "<<prop<<") is:\n"<<prop->show()<<endl;
  //set up the mcmc sampler (assuming mcmc)
  mcmc.setup(Ninit,*like,prior,*prop,output_precision);


  //Prepare for chain output
  ss<<"gle_"<<outname;
  string base=ss.str();

  //Loop over Nchains
  for(int ic=0;ic<Nchain;ic++){
    bayes_sampler *s=s0->clone();
    s->initialize();
    s->run(base,ic);
    {//For exact backward compatibility we need to override the command-line flag -poly and always set integrate=true for the analysis
      //Here we try first working exclusively with the new code ML_photometry_likelihood, rather than MLFitProb...
      //As coded here this is not very general...
      GLensBinary alens;
      alens.Optioned::addOptions(opt,"");
      alens.setup();
      alens.set_integrate(true);//This line is the point of remaing all this.
      cout<<"alens:"<<alens.print_info()<<endl;
      ML_photometry_signal asignal(traj, &alens);
      asignal.Optioned::addOptions(opt,"");
      asignal.setup();
      //asignal.set_tstartHACK(tstart);
      ML_photometry_likelihood alike(&space, data, &asignal, &prior);
      cout<<"alike="<<&alike<<endl;
      alike.Optioned::addOptions(opt,"");
      alike.setup();
      alike.defWorkingStateSpace(space);
      s->analyze(base,ic,Nsigma,Nbest,alike);
    }
    delete s;
  }
  
  //Dump summary info
  cout<<"best_post "<<like->bestPost()<<", state="<<like->bestState().get_string()<<endl;
  }

//An analysis function defined below.
void dump_view(const string &outname, bayes_data &data, ML_photometry_signal &signal, bayes_likelihood &like,state &s,double tstart,double tend,int nsamples){
  //The report includes:
  // 1. The lens magnification map
  // 2. The observer trajectory curve
  // 3. The model lightcurve
  // (anything else?)
  ostringstream ss;
  
  //magnification map
  ss.str("");ss<<outname<<"_mmap.dat";
  dump_mag_map(ss.str(), data, signal, s, tstart, tend, nsamples);

  //trajectory
  vector<double>times;
  times=data.getLabels();
  double tref=data.getFocusLabel(true);

  //data grid trajectory
  ss.str("");ss<<outname<<"_d_traj.dat";
  //dump_trajectory(ss.str(),odata,s,0,0);  
  {
    ofstream out(ss.str());
    signal.dump_trajectory(out, s, times, tref);
  }
  //fine trajectory
    ss.str("");ss<<outname<<"_traj.dat";
  {
    double delta_t=(tend-tstart)/(nsamples-1.0);
    for(int i=0;i<nsamples;i++){
      double t=tstart+i*delta_t;
      times.push_back(t);
    }
    //dump_trajectory(ss.str(),odata,s,tstart,tend,nsamples*4);//use 4 times mm_samples   
    ofstream out(ss.str());
    signal.dump_trajectory(out, s, times, tref);
  }

  //light curve
  ss.str("");ss<<outname<<"_lcrv.dat";
  dump_lightcurve(ss.str(),like,s,tstart,tend,nsamples*4);//use 4 times mm_samples     

  //data light curve
  ss.str("");ss<<outname<<"_d_lcrv.dat";
  dump_lightcurve(ss.str(),like,s,0,0);  

};

void dump_mag_map(const string &outname, bayes_data &data,ML_photometry_signal &signal, state &s,double tstart,double tend,int nsamples,int cent,bool output_nlens){
  ofstream out(outname);
  if(tend<=tstart)data.getDomainLimits(tstart,tend);
  Point LLp(0,0), URp(0,0);
  signal.getWindow(s, LLp, URp, tstart, tend, cent);  
  cout<<"dump mag map: LL=("<<LLp.x<<","<<LLp.y<<") UR=("<<URp.x<<","<<URp.y<<")"<<endl;
  GLens *lens=signal.clone_lens(s);
  cout<<"lens="<<lens->print_info();
  lens->writeMagMap(out, LLp, URp, nsamples, output_nlens);
  delete lens;
};

///Dump the lightcurve
void dump_lightcurve(const string &outname,bayes_likelihood &like,state &s,double tstart,double tend,int nsamples){
  ofstream out(outname);
  if(tend-tstart>0){
    like.writeFine(out,s,nsamples,tstart,tend);
  } else {
    like.write(out,s);
  }
  out<<endl;
};

