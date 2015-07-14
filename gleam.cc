//Gravitational lens event analysis machine
//Written by John G Baker NASA-GSFC (2014)

#include "mlfit.hh"
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
bool integrate;
int Nevery;
const double MaxAdditiveNoiseMag=22;
int output_precision;
double mm_lens_rWB;

//Analysis functions defined below.
void dump_view(const string &outname,MLdata&data,state &s,double tstart,double tend, int nsamples=301);
void dump_mag_map(const string &outname,MLdata&data,state &s,double tstart,double tend, int nsamples=301, int cent=-2,bool output_nlens=false);
void dump_trajectory(const string &outname,MLdata&data,state &s,double tstart,double tend,int nsamples=301);
void dump_lightcurve(const string &outname,MLdata&data,state &s,double tstart,double tend,int nsamples=301);


///Our likelihood function class, defined below
///Perhaps move this into mlfit?

class MLFitProb: public bayes_likelihood{
  MLdata * odata;
  sampleable_probability_function * prior;
  int count;
  double total_eval_time;
  stateSpaceTransformND noise_trans{2,{"I0","Fn"},{"I0","Mn"},[](vector<double>&v){return vector<double>({v[0],v[0]-2.5*log10(v[1])});}};
  bool do_additive_noise;
 public:
  double best_post;
  state best;

  virtual ~MLFitProb(){};//Avoid GCC warnings
  MLFitProb(stateSpace *sp, MLdata *data, bayes_data *d,bayes_signal *s,sampleable_probability_function *prior=nullptr):odata(data),prior(prior),bayes_likelihood(sp,d,s){
    best=state(sp,Npar);
    reset();
    //void setup(){    
    //if(optSet("additive_noise"))useAdditiveNoise();
    set_like0();
    //cout<<"setup:options="<<reportOptions()<<endl;
    //cout<<"setup:do_additive_noise="<<(do_additive_noise?"true":"false")<<endl;
  };
  void reset(){
    best_post=-INFINITY;
    best=best.scalar_mult(0);
    count=0;
    total_eval_time=0;
  }
  state bestState(){return best;};
  double bestPost(){return best_post;};

  double evaluate_log(state &s){
    //#pragma omp critical
    //cout<<"Evaluating likelihood for params"<<s.get_string()<<endl;
    valarray<double>params=s.get_params();
    //clock_t tstart=clock();
    double tstart=omp_get_wtime();
    //double result=odata->evaluate_log(params,integrate);
    //cout<<"loglike="<<result<<"<="<<best_post<<endl;   
    //result=log_chi_squared(s);
    double result=log_chi_squared(s);
    double post=result;
    if(prior)post+=prior->evaluate_log(s);
    //cout<<"loglike="<<result<<"<="<<best_post<<endl;   
    //clock_t tend=clock();
    //double eval_time = (tend-tstart)/(double)CLOCKS_PER_SEC;
    double tend=omp_get_wtime();
    double eval_time = tend-tstart;
    #pragma omp critical
    {     
      total_eval_time+=eval_time;
      count++;
      if(0==count%Nevery)
	cout<<"eval_time = "<<eval_time<<"  result="<<result<<" mean = "<<total_eval_time/count<<endl; 
      if(post>best_post){
        best_post=post;
        best=state(s);
      }
      if(!isfinite(result)){
        cout<<"Whoa dude, loglike is NAN! What's up with that?"<<endl;
        cout<<"params="<<s.get_string()<<endl;
	result=-INFINITY;
      }
    }
    return result;
  };
  //string print_info(){
  //ostringstream s;
  //s<<"MLFitProb:data["<<i<<"]:\n"<<data->print_info();
  //return s.str();
  //};
  void defWorkingStateSpace(const stateSpace &sp){
    if(optSet("additive_noise"))do_additive_noise=true;
    else do_additive_noise=false;
    haveSetup();
    checkSetup();//Call this assert whenever we need options to have been processed.
    haveWorkingStateSpace();
    checkPointers();
    signal->defWorkingStateSpace(sp);
    //backward compatible hack    
    if(!do_additive_noise){
      stateSpace st=noise_trans.transform(sp);
      //cout<<"mlfitprob:defWSS: transformed space is:\n"<<sp.show()<<endl;
      data->defWorkingStateSpace(st);
    } else {
      data->defWorkingStateSpace(sp);
    }    
  };
  state transformDataState(const state &s)const{
    if(!do_additive_noise){
      state st=noise_trans.transformState(s);
      //cout<<"Transforming data state from:"<<s.show()<<"\nto:"<<st.show()<<endl;
      return st;
    }
    return s;
  };
  state transformSignalState(const state &s)const{return s;};

  void addOptions(Options &opt,string s=""){Optioned::addOptions(opt,s);};
  
  stateSpace getObjectStateSpace()const{return stateSpace();};
  void write(ostream &out,state &st){
    debug_signal=true;
    odata->write(out,st);
    //oldwrite(out,st);
  };
  void writeFine(ostream &out,state &st){
    double nsamples=0,tstart=0,tend=0;
    debug_signal=false;
    getFineGrid(nsamples,tstart,tend);
    odata->write(out,st,nsamples,tstart,tend);
    //oldwrite(out,st,nsamples,tstart,tend);
  };
  void getFineGrid(double & nfine, double &tfinestart, double &tfineend)const{
    //nfine=odata->size()*2;
    nfine=data->size()*2;
    double t0,twidth;
    double tstart,tend;
    //odata->getTimeLimits(tstart,tend);
    data->getDomainLimits(tstart,tend);
    //t0=odata->getPeakTime();
    t0=data->getFocusLabel();
    double finewidth=1.5;
    tfinestart=t0-(-tstart+tend)*finewidth/2.0;
    tfineend=t0+(-tstart+tend)*finewidth/2.0; //tfine range is twice data range centered on t0
    twidth=10;//look mostly within 10 days of peak;
    cout<<"tfs="<<tfinestart<<" < ts="<<tstart<<" < t0="<<t0<<" < te="<<tend<<" < tfe="<<tfineend<<endl;
    cout<<"nfine="<<nfine<<endl;
  };

  void oldwrite(ostream &out, state&st, int nsamples=-1, double tstart=0, double tend=0){
    checkPointers();
    vector<double>times;
    if(nsamples<0)
      times=data->getLabels();
    else {
      double delta_t=(tend-tstart)/(nsamples-1);
      for(int i=0;i<nsamples;i++){
	double t=tstart+i*delta_t;
	times.push_back(t);
      }
    }
    //double tpk=getPeakTime();
    double tpk=data->getFocusLabel();
    double time0=data->getFocusLabel(true);
    
    //backward compatibility hack
    const int idx_Fn=2,idx_I0=0;
    double noise_lev=st.get_param(idx_Fn);
    double I0=st.get_param(idx_I0);
    double noise_mag=I0-2.5*log10(noise_lev);
    if(do_additive_noise)noise_mag=noise_lev;
    cout<<"I0,noise_lev,integrate="<<I0<<","<<noise_lev<<","<<integrate<<endl;
    cout<<"nsamples,tstart,tend="<<nsamples<<" "<<tstart<<" "<<tend<<endl;
    //endhack
    
    I0=st.get_param(idx_I0);
    vector<double> model=signal->get_model_signal(st,times);
    vector<double> dmags=data->getDeltaValues();
    vector<double> dvar=getVariances(st);

    if(nsamples<0){
      for(int i=0;i<times.size();i++){
	double S=dvar[i];
	//if(i<10)cout<<"i="<<i<<"  S="<<S<<endl;
	if(i==0)
	  out<<"#t"<<" "<<"t_vs_pk" 
	     <<" "<<"data_mag"<<" "<<"model_mag"
	     <<" "<<"data_err"<<" "<<"model_err"
	     <<endl;
	out<<times[i]+time0<<" "<<times[i]-tpk
	   <<" "<<data->getValue(i)<<" "<<model[i]
	   <<" "<<dmags[i]<<" "<<sqrt(S)
	   <<endl;
      }
    } else {
      for(int i=0;i<times.size();i++){
	//if(i<10)cout<<"i="<<i<<"  S="<<S<<endl;
	double rtS=pow(10.0,0.4*(-noise_mag+model[i]));
	double t=times[i];
	if(i==0)
	  out<<"#t"<<" "<<"t_vs_pk" 
	     <<" "<<"model_mag"<<" "<<"model_extra_err"
	     <<endl;
	out<<t+time0<<" "<<t-tpk
	   <<" "<<model[i]<<" "<<rtS
	   <<endl;
      }
    }
  };


};

int ref_dir;//use to parameter for marginalization.
//Function to estimate the marginalized parameter value
double for_peak_position(state s){
  vector<double>pars=s.get_params_vector();
  return pars[ref_dir];
};

//***************************************************************************************8
//main test program
int main(int argc, char*argv[]){
  string datafile;
  const int NparRead=Npar; 

  //Create the sampler
  ptmcmc_sampler mcmc;
  bayes_sampler *s0=&mcmc;
  //Create the model components
  GLensBinary binarylens;
  Trajectory linear_trajectory(Point(0,0), Point(1,0));
  Trajectory *traj=&linear_trajectory;
  GLens *lens=&binarylens;
  ML_OGLEdata data;
  ML_photometry_signal signal(traj, lens);

  Options opt;

  s0->addOptions(opt,"");
  lens->addOptions(opt,"");
  data.addOptions(opt,"");
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
  opt.add(Option("mm_samples","Number of samples in magmap (default 300)","300"));
  opt.add(Option("mm_nimage","Include number of images magmap"));
  opt.add(Option("precision","Set output precision digits. (Default 13).","13"));
  opt.add(Option("mm_lens_rWB","In magmap mode, radial scale outside which we use the WideBinary lens inversion (neg for none). Use for quality control. (Default 5).","5"));

  int Nlead_args=2;
  
  bool parseBAD=opt.parse(argc,argv);
  if(parseBAD||(argc != Nlead_args+1 && argc!=Nlead_args+1+Npar && argc!=5)) {
    cout << "You gave " << argc-1 << " arguments" << endl;
    cout << "Usage:\n gleam [-options=vals] data_file_name output_name [ I0 Fs Fn logq logL r0 phi tE tpass ]" << endl;
    cout << "Or:\n gleam -magmap [-options=vals] output_name logq logL width" << endl;
    cout <<opt.print_usage()<<endl;
    return 1;
  }
    
  //options
  //DEV: take options from sub-elements rather than filling them all here. eg call "add_options"
  //DEV: Subelements:
  //DEV:   -Model (Now just one option)
  //DEV:   -Bayes-estimator (Now just MCMC)
  //DEV:   -Data?
  cout<<"flags=\n"<<opt.report()<<endl;
  
  double nburn_frac,Tmax,Fn_max,tE_max,tcut,seed;
  int Nchain;
  int mm_center,mm_samples,save_every;
  bool do_magmap,view;
  int Nsigma=1;
  int Nbest=10;
  *s0->optValue("nevery")>>Nevery;
  view=opt.set("view");
  istringstream(opt.value("nchains"))>>Nchain;
  istringstream(opt.value("seed"))>>seed;
  //if seed<0 set seed from clock
  if(seed<0)seed=time(NULL);
  istringstream(opt.value("tcut"))>>tcut;
  do_magmap=opt.set("magmap");
  istringstream(opt.value("mm_center"))>>mm_center;
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
    datafile="";
    outname=argv[1];
    if(Npar<3){
      cout<<"We expected Npar>=3"<<endl;
      exit(1);
    }
    for(int i=0;i<3;i++)stringstream(argv[i+2])>>params[i];  
    have_pars0=true;
  } else if(argc>=1+Nlead_args){
    datafile=argv[1];
    outname=argv[2];
    if(argc>1+Nlead_args){
      for(int i=0;i<NparRead;i++)stringstream(argv[i+1+Nlead_args])>>params[i];
      have_pars0=true;
    }
  }

  //report
  cout.precision(output_precision);
  cout<<"\noutname = '"<<outname<<"'"<<endl;
  cout<<"seed="<<seed<<endl; 
  cout<<"integrate="<<(integrate?"true":"false")<<endl;
  cout<<"Running on "<<omp_get_max_threads()<<" thread"<<(omp_get_max_threads()>1?"s":"")<<"."<<endl;

  ProbabilityDist::setSeed(seed);
  globalRNG.reset(ProbabilityDist::getPRNG());//just for safety to keep us from deleting main RNG in debugging.
  //cout<<"globalRNG="<<globalRNG.get()<<endl;
  
  //Set up ML objects
  data.setup(datafile);
  signal.setup();


  //Create the data object (This block is setting up the model.  If/when we separate model from data, then data obj may be defined below)
  cout<<"OGLE data file='"<<datafile<<"'"<<endl;
  OGLEdata odata(datafile);
  double r0s=6.0,q0;
  istringstream(opt.value("q0"))>>q0;
  bool use_remapped_r0,use_remapped_q,use_log_tE,use_additive_noise=false;
  use_remapped_r0=opt.set("remap_r0");
  use_remapped_q=opt.set("remap_q");
  use_additive_noise=opt.set("additive_noise");
  use_log_tE=opt.set("log_tE");
  if(use_remapped_r0){
    odata.remap_r0(2.0);
    r0s=1.0;
  }
  if(use_remapped_q){
    odata.remap_q(q0);
  }
  if(use_additive_noise){
    odata.useAdditiveNoise();    
  }
  if(use_log_tE){
    odata.use_log_tE();
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
    double qpar,Lpar,widthpar;
    qpar=params[0];
    Lpar=params[1];
    widthpar=params[2];
    params[3]=qpar;
    params[4]=Lpar;
    params[5]=widthpar;
    params[7]=1;
    //other params are undefined and should be unused...
    state s(&space,params);
    //magnification map
    ss.str("");ss<<outname<<"_mmap.dat";
    dump_mag_map(ss.str(),odata,s,0,0,mm_samples,mm_center);    
    if(opt.set("mm_nimage")){
      //debugint=true;
      ss.str("");ss<<outname<<"_nmap.dat";    
      dump_mag_map(ss.str(),odata,s,0,0,mm_samples,mm_center,true);      
    }
    exit(0);
  }    

  //Prune data
  odata.cropBefore(tcut);
  cout<<"Ndata="<<odata.size()<<endl;
  double t0,twidth;
  //define time range:
  double tstart,tend;
  odata.getTimeLimits(tstart,tend);
  t0=odata.getPeakTime();
  double finewidth=1.5;
  double tfinestart=t0-(-tstart+tend)*finewidth/2.0;
  double tfineend=t0+(-tstart+tend)*finewidth/2.0; //tfine range is twice data range centered on t0
  twidth=10;//look mostly within 10 days of peak;
  cout<<"tfs="<<tfinestart<<" < ts="<<tstart<<" < t0="<<t0<<" < te="<<tend<<" < tfe="<<tfineend<<endl;

  //Initial parameter state
  state instate(&space,params);
  if(have_pars0){
    cout<<"Input parameters:"<<instate.show()<<endl;
  }

  ///At this point we are ready for analysis in the case that we are asked to view a model
  ///Note that we still have needed the data file to create the OGLEdata object, and concretely
  ///to set the domain.  This could be changed...
  if(view){
    if(have_pars0){
      cout<<"Producing report on the model with specified parameters."<<endl;
      dump_view(outname,odata,instate,tfinestart,tfineend,mm_samples);
    } else {
      cout<<" The -view option requires that parameters are provided."<<endl;
    }
  }

  //Set the prior:
  //Eventually this should move to the relevant constitutent code elements where the params are given meaning.
  const int uni=mixed_dist_product::uniform, gauss=mixed_dist_product::gaussian, pol=mixed_dist_product::polar; 
  //                                     I0      Fs     Fn          logq*   logL       r0     phi      tE     tpass  
  valarray<double>    centers((dlist){ 18.0,    0.5,   0.5*Fn_max,   0.0,   0.0,  r0s/2.0,   M_PI, tE_max/2,  t0     });
  valarray<double> halfwidths((dlist){  5.0,    0.5,   0.5*Fn_max,   2.0,   2.0,  r0s/2.0,   M_PI, tE_max/2,  twidth });
  valarray<int>         types((ilist){gauss,    uni,   uni,        gauss, gauss,      uni,    uni,      uni,  gauss  });
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
  mixed_dist_product prior(&space,types,centers,halfwidths);
  cout<<"Prior is:\n"<<prior.show()<<endl;

  //Set the likelihood
  //FIXME some of this should move up before addOptions 
  lens->setup();
  bayes_likelihood *llike=nullptr;
  bayes_likelihood *like=nullptr;
  llike = new MLFitProb(&space,&odata,&data,&signal,&prior);
  llike->addOptions(opt,"");
  ML_photometry_likelihood mpl(&space, &data, &signal, &prior);
  mpl.addOptions(opt,"");
  mpl.setup();
  like = &mpl;
  llike->defWorkingStateSpace(space);
  like->defWorkingStateSpace(space);
  ///At this point we are ready for analysis in the case that we are asked to view a model
  ///Note that we still have needed the data file to create the OGLEdata object, and concretely
  ///to set the domain.  This could be changed...
  if(view){
    if(have_pars0){
    double nsamples,tfinestart,tfineend;
    llike->getFineGrid(nsamples,tfinestart,tfineend);
      cout<<"Producing report on the model with specified parameters."<<endl;
      dump_view(outname,odata,instate,tfinestart,tfineend,mm_samples);
    } else {
      cout<<" The -view option requires that parameters are provided."<<endl;
    }
  }

  //Just a little reporting
  if(have_pars0){
    double ll=llike->evaluate_log(instate);
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
  mcmc.setup(Ninit,*llike,prior,*prop,output_precision);


  //Prepare for chain output
  ss<<"gle_"<<outname;
  string base=ss.str();

  //Loop over Nchains
  for(int ic=0;ic<Nchain;ic++){
    bayes_sampler *s=s0->clone();
    s->initialize();
    s->run(base,ic);
    //s->analyze(base,ic,Nsigma,Nbest,*llike);
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
      ML_photometry_likelihood alike(&space, &data, &asignal, &prior);
      cout<<"alike="<<&alike<<endl;
      alike.Optioned::addOptions(opt,"");
      alike.setup();
      alike.defWorkingStateSpace(space);
      s->analyze(base,ic,Nsigma,Nbest,alike);
    }
  }
  
  //Dump summary info
  cout<<"best_post "<<llike->bestPost()<<", state="<<llike->bestState().get_string()<<endl;
}

//An analysis function defined below.
void dump_view(const string &outname,MLdata&data,state &s,double tstart,double tend,int nsamples){
  //The report includes:
  // 1. The lens magnification map
  // 2. The observer trajectory curve
  // 3. The model lightcurve
  // (anything else?)
  ostringstream ss;
  
  //magnification map
  ss.str("");ss<<outname<<"_mmap.dat";
  dump_mag_map(ss.str(),data,s,tstart,tend,nsamples);  

  //trajectory
  ss.str("");ss<<outname<<"_traj.dat";
  dump_trajectory(ss.str(),data,s,tstart,tend,nsamples*4);//use 4 times mm_samples   

  //data trajectory
  ss.str("");ss<<outname<<"_d_traj.dat";
  dump_trajectory(ss.str(),data,s,0,0);  

  //light curve
  ss.str("");ss<<outname<<"_lcrv.dat";
  dump_lightcurve(ss.str(),data,s,tstart,tend,nsamples*4);//use 4 times mm_samples     

  //data light curve
  ss.str("");ss<<outname<<"_d_lcrv.dat";
  dump_lightcurve(ss.str(),data,s,0,0);  

};

//This is like glensscan.cc
void dump_mag_map(const string &outname,MLdata&data,state &s,double tstart,double tend,int nsamples,int cent,bool output_nlens){
  debug=false;
  valarray<double>params=s.get_params();  
  bool use_r0_width=false;
  double I0,Fs,noise_frac,q,L,r0,phi,tE,tmax;
  data.get_model_params(params, I0,Fs,noise_frac,q,L,r0,phi,tE,tmax);
  GLensBinary lens(q,L);
  lens.set_WideBinaryR(mm_lens_rWB);
  double x0=0,y0=0,width=0,wx,wy,xcm  =  (q/(1.0+q)-0.5)*L;

  //Get square bounding box size
  cout<<"q,L,cent="<<q<<", "<<L<<", "<<cent<<endl;
  if(cent>-2){//signal to use r0 as map width and center on {rminus,CoM,rplus}, when cent={0,1,2}
    //We have to hand code the CoM info
    width=r0;
    y0=-width/2;
    switch(cent){
    case -1:
      x0=-0.5*L;
    case 0:
      x0=xcm;
    case 1:
      x0=0.5*L;
    }
  } else {
    double tleft=(tstart-tmax)/tE,tright=(tend-tmax)/tE;
    if(tend>tstart){
      tleft=(tstart-tmax)/tE;
      tright=(tend-tmax)/tE;
    } else { //use data range
      data.getTimeLimits(tleft,tright);
      tleft=(tleft-tmax)/tE;    
      tright=(tright-tmax)/tE;    
    }  
    Trajectory traj(data.get_trajectory(q,L,r0,phi,tleft));
    Point pstart=traj.get_obs_pos(tleft);
    Point pend=traj.get_obs_pos(tright);
    cout<<"making mag-map between that fits points: ("<<pstart.x<<","<<pstart.y<<") and ("<<pend.x<<","<<pend.y<<")"<<endl;
    width=wx=abs(pstart.x-pend.x);
    wy=abs(pstart.y-pend.y);
    if(wy>width)width=wy;
    x0=pstart.x;
    if(pend.x<x0)x0=pend.x;
    y0=pstart.y;
    if(pend.y<y0)y0=pend.y;
    width+=1.0;//margin
    y0=y0-(width-wy)/2.0;
    x0=x0-(width-wx)/2.0;
  }
  cout<<"x0,y0,width="<<x0<<", "<<y0<<", "<<width<<endl;
    
  double dx=width/(nsamples-1);    
  cout<<"mag-map ranges from: ("<<x0<<","<<y0<<") to ("<<x0+width<<","<<y0+width<<") stepping by: "<<dx<<endl;
  ofstream out(outname);
  double ten2prec=pow(10,output_precision-2);
  out.precision(output_precision);
  out<<"#x-xcm  y  magnification"
    //<<"  Nimages"
     <<endl;
  for(double y=y0;y<=y0+width;y+=dx){
    //cout<<"y="<<y<<endl;
    Trajectory traj(Point(x0,y), Point(1,0), width, dx);
    vector<int> indices;
    vector<double> times,mags;
    vector<vector<Point> >thetas;
    //debugint=true;
    lens.compute_trajectory(traj,times,thetas,indices,mags,integrate);
    //debugint=false;
    int ilast=-1;
    for(int i : indices){
      Point b=traj.get_obs_pos(times[i]);
      double mtruc=floor(mags[i]*ten2prec)/ten2prec;
      out.precision(output_precision);
      out<<b.x-xcm<<" "<<b.y<<" "<<setiosflags(ios::scientific)<<mtruc<<setiosflags(ios::fixed);
      if(output_nlens)out<<" "<<thetas[i].size();
      out<<endl;
      ilast=i;
    }
    out<<endl;
  }	  
};

///Dump the trajectory
///And additional information about the lens evaluation.
void dump_trajectory(const string &outname,MLdata&data,state &s,double tstart,double tend,int nsamples){
  valarray<double>params=s.get_params();  
  double I0,Fs,noise_frac,q,L,r0,phi,tE,tmax;
  data.get_model_params(params, I0,Fs,noise_frac,q,L,r0,phi,tE,tmax);
  double dt=(-tstart+tend)/tE/nsamples;
  double tleft=(tstart-tmax)/tE,tright=(tend-tmax)/tE;
  double tpk=data.getPeakTime();
  double time0=data.getPeakTime(true);
  //Get square bounding box size
  cout<<"Model params: I0  Fs  noise_frac; q  L  r0  phi tE tmax\n"
      <<I0<<" "<<Fs<<" "<<noise_frac<<"; "<<q<<" "<<L<<" "<<r0<<" "<<phi<<" "<<tE<<" "<<tmax<<endl;
  cout<<"making traj("<<q<<" "<<L<<" "<<r0<<" "<<phi<<")"<<endl;
  cout<<" vx="<<cos(phi)<<", vy="<<sin(phi)<<endl;
  cout<<" p0x="<<-r0*sin(phi)<<", p0y="<<r0*cos(phi)<<endl;
  cout<<" xcm  =  "<<(q/(1.0+q)-0.5)*L<<endl;
  double nu=1/(1+q),nu_pos=1-nu;

  vector<double>times;
  vector<double>xtimes,modelmags; 
  vector<int> indices;
  //cout<<"tleft="<<tleft<<" tright="<<tright<<" dt="<<dt<<endl;
  if(dt>0){//if no time range use data times...
    cout<<"traj: Setting times provided: from "<<tleft<<" to "<<tright<<endl;
    for(double t=tleft;t<=tright;t+=dt){
      times.push_back(t);
    }
  } else {
    cout<<"traj: Setting times from data"<<endl;
    times=data.getTimes();
    for(double &t:times)t=(t-tmax)/tE;
  }

  Trajectory traj(data.get_trajectory(q,L,r0,phi,times[0]));
  traj.set_times(times,0);
  cout<<"times range from "<<traj.t_start()<<" to "<<traj.t_end()<<endl;

  ofstream out(outname);
  GLensBinary lens(q,L);
  vector<vector<Point>> allthetas;
  
  //debug=true;
  lens.compute_trajectory(traj,xtimes,allthetas,indices,modelmags,true);
  //debug=false;

  out.precision(output_precision);
  cout<<traj.print_info()<<endl;
  //out<<traj.print_info()<<endl;
  out<<"#"<<s.get_string()<<endl;
  out<<"#nu_pos="<<nu_pos<<" nu_neg="<<nu<<" L="<<L<<" r0="<<r0<<" phi="<<phi<<endl;
  out<<"#1.t   2. t_rel  3.x   4.y   5.magnif   6.n_img   7.magnif(wide)   8.n_img   9.magnif(int)   10. n_img 11.rpos  12.rneg 13-22.imgsXY 23-28.wide-imgsXY 29-38.int-imgsXY"<<endl;
  //cout<<tstart<<"< t <"<<tend<<" dt="<<dt*tE<<" tE="<<tE<<endl;
  //cout<<tleft<<"< t* <="<<tright<<endl;
  //cout<<"xtimes.size()="<<xtimes.size()<<endl;
  int i=0;
  //for(double t=tleft;t<=tright;t+=dt){
  for(int ii:indices){
    //cout<<"indices["<<i<<"]="<<ii<<endl;
    double t=xtimes[ii];
    //cout<<"t="<<t<<endl;
    Point p=traj.get_obs_pos(t);
    vector<Point> thetas=lens.invmap(p);
    vector<Point> wthetas=lens.invmapWideBinary(p);
    double x1=p.x-L/2,x2=p.x+L/2,r1sq=x1*x1+p.y*p.y,r2sq=x2*x2+p.y*p.y;
    double cpos=nu_pos/sqrt(r1sq),cneg=nu/sqrt(r2sq),xcm  =  (q/(1.0+q)-0.5)*L;
    
    //out<<t<<" "<<xtimes[indices[i]]<<" "<<p.x-xcm<<" "<<p.y<<" "<<lens.mag(thetas)<<" "<<thetas.size()<<" "<<lens.mag(wthetas)<<" "<<wthetas.size()
    // <<" "<<modelmags[indices[i]]<<" "<<allthetas[indices[i]].size();
    double Ival = I0 - 2.5*log10(Fs*modelmags[indices[i]]+1-Fs);
    out<<t+time0<<" "<<t-tpk<<" "<<p.x-xcm<<" "<<p.y<<" "<<lens.mag(thetas)<<" "<<thetas.size()<<" "<<lens.mag(wthetas)<<" "<<wthetas.size()
       <<" "<<modelmags[ii]<<" "<<allthetas[ii].size();
    //out<<t<<" "<<xtimes[indices[i]]<<" "<<p.x<<" "<<p.y<<" "<<lens.mag(wthetas)<<" "<<wthetas.size();

    for(int j=0;j<5;j++){
      if(j<thetas.size())out<<" "<<thetas[j].x<<" "<<thetas[j].y;
      else out<<" 0 0";
      }
      
    for(int j=0;j<3;j++){
      if(j<wthetas.size())out<<" "<<wthetas[j].x<<" "<<wthetas[j].y;
      //{
	//Point betaj=lens.map(wthetas[j]);
	//out<<" "<<wthetas[j].x<<"->"<<betaj.x-p.x
	// <<" "<<wthetas[j].y<<"->"<<betaj.y-p.y;
      //}
      else out<<" 0 0";
    }
    
    for(int j=0;j<5;j++){
      //if(j<allthetas[indices[i]].size())out<<" "<<allthetas[indices[i]][j].x<<" "<<allthetas[indices[i]][j].y;
      if(j<allthetas[ii].size())out<<" "<<allthetas[ii][j].x<<" "<<allthetas[ii][j].y;
      else out<<" 0 0";
      }
    out<<endl;
    i++;
  }
};

///Dump the lightcurve
void dump_lightcurve(const string &outname,MLdata&data,state &s,double tstart,double tend,int nsamples){
  ofstream out(outname);
  int nfine;
  if(tend-tstart>0)
    nfine=nsamples;
  else 
    nfine=-1;
  out.precision(output_precision);    
  vector<double> p=s.get_params_vector();
  out<<"#"<<s.get_string()<<endl;
  data.write(out,p,true,nfine,tstart,tend);
  out<<endl;
};

