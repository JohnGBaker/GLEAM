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
const double MaxAdditiveNoiseMag=22;
int output_precision;
double mm_lens_rWB;

//Analysis functions defined below.
void dump_view(const string &outname, bayes_data &data, ML_photometry_signal &signal, bayes_likelihood &like,state &s,double tstart,double tend,int nsamples);
void dump_mag_map(const string &outname,bayes_data &data,ML_photometry_signal &signal, state &s,double tstart,double tend,int nsamples=301);//,int cent=-2,bool output_nlens=false);
void dump_lightcurve(const string &outname,bayes_likelihood&like,state &s,double tstart,double tend,int nsamples=301);


//***************************************************************************************8
//main test program
int main(int argc, char*argv[]){

  //Create the sampler
  ptmcmc_sampler mcmc;
  bayes_sampler *s0=&mcmc;
  //Create the model components
  GLens singlelens;
  GLensBinary binarylens;
  bayes_component_selector lenses(vector<bayes_component*>({&binarylens,&singlelens}));
  Trajectory linear_trajectory;
  Trajectory *traj=&linear_trajectory;
  //GLens *lens=&binarylens;
  GLens *lens;
  bool do_mock=false;
  
  ///So far there are 2 tested types ML_OGLE_data and ML_generic_data while ML_photometry_mock_data is
  ///under development.  [Thinking ahead, we may later have
  ///astrometric data as well. These might be similarly added, but they will involve a different type of
  ///signal model.  We may need to be able to use vector-valued bayes_data::values with a vector of type
  ///channel-labels that can be checked.  Maybe the channels are just strings, but if there is important 
  ///meta-information regarding the channels, then these could become instances of a channel class.]
  Options opt(false);
  ML_photometry_data::addStaticOptions(opt);
  //For the first-pass option check, we first make a copy of argc/argv
  int ac=argc;char* av[argc];
  char * filename=NULL;
  for(int i=0;i<argc;i++){av[i]=argv[i];};
  opt.parse(ac,av,false);
  //select lens
  lens=dynamic_cast<GLens*>(lenses.select(opt));
  //now select the data obj
  ML_photometry_data *data;
  if(opt.set("OGLE_data"))
    data=new ML_OGLEdata();
  else if(opt.set("gen_data"))
    data=new ML_generic_data();
  else if(opt.set("mock_data")){
    data=new ML_mock_data();
    do_mock=true;
  } else {
    //for backward compatibility [deprecated] default is to assume OGLE data and try to read the data from a file named in the (extra) first argument
    if(argc>=1){
      cout<<"Setting filename from first argument for backward compatibility [deprecated]."<<endl;
      filename=av[1];
      data=new ML_OGLEdata;
    } else { //go on assuming mock data; Except for backward compatibility, this would be the default.
      cout<<"No data file indicated!"<<endl;
      data=new ML_mock_data();
      do_mock=true;
    }
  }
      
  //Eventually want to handle signal polymorphism similarly
  ML_photometry_signal signal(traj, lens);
  bayes_likelihood *like=nullptr;
  ML_photometry_likelihood mpl(data, &signal);
  mpl.addOptions(opt);
  like=&mpl;
  s0->addOptions(opt);
  lens->addOptions(opt);
  traj->addOptions(opt);
  data->addOptions(opt);
  signal.addOptions(opt);

  opt.add(Option("nchains","Number of consequtive chain runs. Default 1","1"));
  opt.add(Option("seed","Pseudo random number grenerator seed in [0,1). (Default=-1, use clock to seed.)","-1"));
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
  if(parseBAD) {
    cout << "Usage:\n gleam [-options=vals] data_file_name output_name [ params ]" << endl;
    cout << "Or:\n gleam -magmap [-options=vals] output_name logq logL width" << endl;
    cout <<opt.print_usage()<<endl;
    return 1;
  }
    
  cout<<"flags=\n"<<opt.report()<<endl;

  //Post parse setup
  lens->setup();  
  traj->setup();  

  //double nburn_frac,Tmax,Fn_max,tE_max,tcut,seed;
  double nburn_frac,Tmax,Fn_max,tcut,seed;
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
  if(seed<0)seed=fmod(time(NULL)/3.0e7,1);
  istringstream(opt.value("tcut"))>>tcut;
  do_magmap=opt.set("magmap");
  istringstream(opt.value("mm_center"))>>mm_center;
  istringstream(opt.value("mm_d0x"))>>mm_d0x;
  istringstream(opt.value("mm_d0y"))>>mm_d0y;
  istringstream(opt.value("mm_samples"))>>mm_samples;
  istringstream(opt.value("precision"))>>output_precision;
  istringstream(opt.value("mm_lens_rWB"))>>mm_lens_rWB;

  //read non-parameter args
  string outname;
  ostringstream ss("");
  if(argc<Nlead_args+1) {
    cout << "You gave " << argc-1 << " arguments. Expecting "<<Nlead_args<< endl;
    cout << "Usage:\n gleam [-options=vals] data_file_name output_name [ params ]" << endl;
    cout << "Or:\n gleam -magmap [-options=vals] output_name logq logL width" << endl;
    cout <<opt.print_usage()<<endl;
    return 1;
  }
  outname=argv[1];

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
  like->setup();
  cout<<"Ndata="<<data->size()<<endl;

  //Get the space/prior for use here
  stateSpace space;
  shared_ptr<const sampleable_probability_function> prior;  
  space=*like->getObjectStateSpace();
  cout<<"like.nativeSpace=\n"<<space.show()<<endl;
  prior=like->getObjectPrior();
  cout<<"Prior is:\n"<<prior->show()<<endl;
  valarray<double> scales;prior->getScales(scales);

  //Read Params
  int Npar=space.size();
  cout<<"Npar="<<Npar<<endl;
  valarray<double>params;
  bool have_pars0=false;
  if(do_magmap){//different parameters in this case
    Npar=3;
    if(argc!=5) {
      cout << "You gave " << argc-1 << " arguments. Expecting 4."<< endl;
      cout << "Usage:\n gleam [-options=vals] data_file_name output_name [ params ]" << endl;
      cout << "Or:\n gleam -magmap [-options=vals] output_name logq logL width" << endl;
      cout <<opt.print_usage()<<endl;
      return 1;
    }
    params.resize(3);
    for(int i=0;i<3;i++)stringstream(argv[i+2])>>params[i];  
    have_pars0=true;
  } else if(argc==Nlead_args+Npar+1){//Read params as defined via the flags
    params.resize(Npar);
    if(argc>Nlead_args+1){
      for(int i=0;i<Npar;i++)stringstream(argv[i+1+Nlead_args])>>params[i];
      have_pars0=true;
    }
  } else if(!argc==Nlead_args+1){
    cout << "You gave " << argc-1 << " arguments. Expecting "<<Nlead_args<<"."<< endl;
      cout << "Usage:\n gleam [-options=vals] data_file_name output_name [ params ]" << endl;
      cout << "Or:\n gleam -magmap [-options=vals] output_name logq logL width" << endl;
      cout <<opt.print_usage()<<endl;
      return 1;
  }
  //Initial parameter state
  state instate(&space,Npar);
  if(have_pars0){
    instate=state(&space,params);
    cout<<"Input parameters:"<<instate.show()<<endl;
  } else { //If no params give just draw something.  Perhaps only relevant for mock_data
    cout<<"Drawing a state from the prior distribution."<<endl;
    instate=prior->drawSample(*globalRNG);
  }

  //Mock data
  if(do_mock){
    like->mock_data(instate);
    //cout<<"outname is:'"<<outname<<"'"<<endl;
    ss.str("");ss<<outname<<"_mock.dat";
    ofstream outmock(ss.str().c_str());
    outmock.precision(output_precision);    
    outmock<<"#mock data\n#"<<instate.get_string()<<endl;
    like->write(outmock,instate);
    outmock<<endl;
    if(!view)exit(0);
  }
    

  //Handle separate magmap (only) option.
  if(do_magmap){//translate parameters
    //Points referenced in this block refer to *lens frame* consider shifting
    double q,L,width;
    q=pow(10.0,params[0]);
    L=pow(10.0,params[1]);
    width=params[2];
    stateSpace lensSpace=*lens->getObjectStateSpace();
    lens->defWorkingStateSpace(lensSpace);
    state lens_state(&lensSpace,valarray<double>({q,L,0}));
    lens->setState(lens_state);
    Point x0=lens->getCenter(mm_center);
    cout<<"cent="<<mm_center<<" = ("<<x0.x+mm_d0x<<","<<x0.y+mm_d0y<<" xcm="<<lens->getCenter().x<<endl;
    Point pstart(x0.x+mm_d0x-width/2,x0.y+mm_d0y-width/2);
    Point pend(x0.x+mm_d0x+width/2,x0.y+mm_d0y+width/2);
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
      lens->verboseWrite();
      lens->writeMagMap(out, pstart, pend, mm_samples);
    }
    exit(0);
  }    

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
    double lp=prior->evaluate_log(instate);
    cout<<"log-Likelihood at input parameters = "<<ll<<endl;
    cout<<"log-posterior at input parameters = "<<ll+lp<<endl;
  }
  if(view)exit(0);
	  
  //assuming mcmc:
  //Set the proposal distribution 
  int Ninit;
  //proposal_distribution *prop=ptmcmc_sampler::new_proposal_distribution(Npar,Ninit,opt,prior,&halfwidths);
  proposal_distribution *prop=ptmcmc_sampler::new_proposal_distribution(Npar,Ninit,opt,prior.get(),&scales);
  cout<<"Proposal distribution is:\n"<<prop->show()<<endl;
  //set up the mcmc sampler (assuming mcmc)
  mcmc.setup(Ninit,*like,*prior,*prop,output_precision);


  //Prepare for chain output
  //ss<<"gle_"<<outname;
  ss<<outname;
  string base=ss.str();

  //Loop over Nchains
  for(int ic=0;ic<Nchain;ic++){
    bayes_sampler *s=s0->clone();
    s->initialize();
    s->run(base,ic);
    s->analyze(base,ic,Nsigma,Nbest,*like);
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

void dump_mag_map(const string &outname, bayes_data &data,ML_photometry_signal &signal, state &s,double tstart,double tend,int nsamples){//,int cent,bool output_nlens){
  //Points in the routine are in *lens frame* //consider shifting to Trajectory frame.
  ofstream out(outname);
  if(tend<=tstart)data.getDomainLimits(tstart,tend);
  Point LLp(0,0), URp(0,0);
  //signal.getWindow(s, LLp, URp, tstart, tend, cent);  
  signal.getWindow(s, LLp, URp, tstart, tend);  
  cout<<"dump mag map: LL=("<<LLp.x<<","<<LLp.y<<") UR=("<<URp.x<<","<<URp.y<<")"<<endl;
  GLens *lens=signal.clone_lens();
  lens->setState(s);
  out.precision(13);
  cout<<"lens="<<lens->print_info();
  //if(output_nlens)lens->verboseWrite();
  lens->writeMagMap(out, LLp, URp, nsamples);
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

