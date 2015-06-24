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
#include "chain.hh"

using namespace std;

typedef initializer_list<double> dlist;
typedef initializer_list<int> ilist;

bool debug = false;
bool debugint = false;

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
//Other
proposal_distribution* new_proposal_distribution(int & Ninit,const Options &opt, sampleable_probability_function * prior=nullptr, const valarray<double>*halfwidths=nullptr);


///Our likelihood function class, defined below
///Perhaps move this into mlfit?

class MLFitProb: public probability_function{
  MLdata * data;
  sampleable_probability_function * prior;
  int count;
  double total_eval_time;
 public:

  double best_post;
  state best;
  virtual ~MLFitProb(){};//Avoid GCC warnings
  MLFitProb(stateSpace *sp, MLdata *data, sampleable_probability_function *prior=nullptr):data(data),prior(prior),probability_function(sp){
    best=state(sp,Npar);
    reset();
  };
  void reset(){
    best_post=-INFINITY;
    best=best.scalar_mult(0);
    count=0;
    total_eval_time=0;
  }
    
  double evaluate_log(state &s){
    //#pragma omp critical
    //cout<<"Evaluating likelihood for params"<<s.get_string()<<endl;
    valarray<double>params=s.get_params();
    //clock_t tstart=clock();
    double tstart=omp_get_wtime();
    double result=data->evaluate_log(params,integrate);
    double post=result;
    if(prior)post+=prior->evaluate_log(s);
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
      //cout<<"loglike="<<result<<"<="<<maxLike<<endl;   
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

  Options opt;

  opt.add(Option("nchains","Number of consequtive chain runs. Default 1","1"));
  opt.add(Option("seed","Pseudo random number grenerator seed in [0,1). (Default=-1, use clock to seed.)","-1"));
  opt.add(Option("nevery","Frequency to dump chain info. Default=1000.","1000"));
  //ptmcmc opts
  opt.add(Option("save_every","Frequency to store chain info. Default=1.","1"));
  opt.add(Option("nsteps","How long to run the chain. Default=5000.","5000"));
  opt.add(Option("nskip","Only dump every nskipth element. Default=10.","10"));
  opt.add(Option("burn_frac","Portion of chain to disregard as burn-in for some calculations. Default=0.5","0.5"));
  opt.add(Option("pt","Do parallel tempering."));
  opt.add(Option("pt_n","Number of parallel tempering chains. Default 20","20"));
  opt.add(Option("pt_swap_rate","Frequency of parallel tempering swap_trials. Default 0.01","0.01"));
  opt.add(Option("pt_Tmax","Max temp of parallel tempering chains. Default 100","100"));
  opt.add(Option("pt_evolve_rate","Rate at which parallel tempering temps should be allowed to evolve. Default none.","0"));
  opt.add(Option("pt_reboot_rate","Max frequency of rebooting poorly performing parallel tempering chains. Default 0","0"));
  opt.add(Option("pt_reboot_every","How often to test for rebooting poorly performing parallel tempering chains. Default 0","0"));
  opt.add(Option("pt_reboot_grace","Grace period protecting infant instances reboot. Default 0","0"));
  opt.add(Option("pt_reboot_cut","Posterior difference cutoff defining poorly performing parallel tempering chains. Default 100","100"));
  opt.add(Option("pt_reboot_thermal","Temperature dependent cutoff term in defining poorly performing parallel tempering chains. Default 0","0"));
  opt.add(Option("pt_reboot_blindly","Do aggressive random rebooting at some level even if no gaps are found. Default 0","0"));
  opt.add(Option("pt_reboot_grad","Let the reboot grace period depend linearly on temp level with given mean. (colder->longer)"));
  opt.add(Option("prop","Proposal type (0-6). Default=4 (DE with Snooker updates w/o prior draws.)","4"));
  opt.add(Option("de_ni","Differential-Evolution number of initialization elements per dimension. Default=10.","10"));
  opt.add(Option("de_eps","Differential-Evolution gaussian scale. Default=1e-4.","1e-4"));
  opt.add(Option("de_reduce_gamma","Differential Evolution reduce gamma parameter by some factor from nominal value. Default=1.","1"));
  opt.add(Option("de_mixing","Differential-Evolution support mixing of parallel chains."));
  opt.add(Option("de_Tmix","Differential-Evolution degree to encourage mixing info from different temps.(default=300)","300"));
  //end ptmcmc opts
  opt.add(Option("poly","Don't use integration method for lens magnification, use only the polynomial method."));
  opt.add(Option("log_tE","Use log10 based variable (and Gaussian prior with 1-sigma range [0:log10(tE_max)] ) for tE parameter rather than direct tE value."));
  opt.add(Option("remap_r0","Use remapped r0 coordinate."));
  opt.add(Option("remap_q","Use remapped mass-ratio coordinate."));
  opt.add(Option("additive_noise","Interpret Fn as magnitude of additive noise. Fn_max is magnitude of maximum noise level (i.e. minimum noise magnitude)"));
  opt.add(Option("q0","Prior max in q (with q>1) with remapped q0. Default=1e5/","1e5"));
  opt.add(Option("Fn_max","Uniform prior max (min for additive) in Fn. Default=1.0 (18.0 additive)/","1"));
  opt.add(Option("tE_max","Uniform prior max in tE. Default=100.0/","100.0"));
  opt.add(Option("tcut","Cut times before tcut (relative to tmax). Default=-1e20/","-1e20"));
  opt.add(Option("view","Don't run any chains, instead take a set of parameters and produce a set of reports about that lens model."));
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
  double nburn_frac,Tmax,q0,Fn_max,tE_max,tcut,seed,swap_rate,pt_reboot_rate,pt_reboot_cut,pt_reboot_thermal,pt_reboot_blindly,pt_evolve_rate;
  int Nchain,Nstep,Nskip,Nptc,mm_center,mm_samples,save_every,pt_reboot_every,pt_reboot_grace;
  bool parallel_tempering,use_remapped_r0,use_remapped_q,use_log_tE,view,use_additive_noise=false,do_magmap,pt_reboot_grad;
  int Nsigma=1;
  int Nbest=10;
  view=opt.set("view");
  istringstream(opt.value("nchains"))>>Nchain;
  istringstream(opt.value("seed"))>>seed;
  //if seed<0 set seed from clock
  if(seed<0)seed=time(NULL);
  istringstream(opt.value("nevery"))>>Nevery;
  istringstream(opt.value("save_every"))>>save_every;
  istringstream(opt.value("nsteps"))>>Nstep;
  istringstream(opt.value("nskip"))>>Nskip;
  istringstream(opt.value("burn_frac"))>>nburn_frac;
  parallel_tempering=opt.set("pt");
  integrate=!opt.set("poly");
  istringstream(opt.value("pt_n"))>>Nptc;
  istringstream(opt.value("pt_evolve_rate"))>>pt_evolve_rate;
  istringstream(opt.value("pt_reboot_rate"))>>pt_reboot_rate;
  istringstream(opt.value("pt_reboot_every"))>>pt_reboot_every;
  istringstream(opt.value("pt_reboot_grace"))>>pt_reboot_grace;
  istringstream(opt.value("pt_reboot_cut"))>>pt_reboot_cut;
  istringstream(opt.value("pt_reboot_thermal"))>>pt_reboot_thermal;
  istringstream(opt.value("pt_reboot_blindly"))>>pt_reboot_blindly;
  pt_reboot_grad=opt.set("pt_reboot_grad");
  istringstream(opt.value("pt_swap_rate"))>>swap_rate;
  istringstream(opt.value("pt_Tmax"))>>Tmax;
  use_remapped_r0=opt.set("remap_r0");
  use_remapped_q=opt.set("remap_q");
  use_log_tE=opt.set("log_tE");
  use_additive_noise=opt.set("additive_noise");
  istringstream(opt.value("q0"))>>q0;
  istringstream(opt.value("Fn_max"))>>Fn_max;
  istringstream(opt.value("tE_max"))>>tE_max;
  istringstream(opt.value("tcut"))>>tcut;
  do_magmap=opt.set("magmap");
  istringstream(opt.value("mm_center"))>>mm_center;
  istringstream(opt.value("mm_samples"))>>mm_samples;
  istringstream(opt.value("precision"))>>output_precision;
  istringstream(opt.value("mm_lens_rWB"))>>mm_lens_rWB;

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
  
  //Create the data object (This block is setting up the model.  If/when we separate model from data, then data obj may be defined below)
  cout<<"OGLE data file='"<<datafile<<"'"<<endl;
  OGLEdata data(datafile);
  double r0s=6.0;
  if(use_remapped_r0){
    data.remap_r0(2.0);
    r0s=1.0;
  }
  if(use_remapped_q){
    data.remap_q(q0);
  }
  if(use_additive_noise){
    data.useAdditiveNoise();    
  }
  if(use_log_tE){
    data.use_log_tE();
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
    dump_mag_map(ss.str(),data,s,0,0,mm_samples,mm_center);    
    if(opt.set("mm_nimage")){
      //debugint=true;
      ss.str("");ss<<outname<<"_nmap.dat";    
      dump_mag_map(ss.str(),data,s,0,0,mm_samples,mm_center,true);      
    }
    exit(0);
  }    

  //Prune data
  data.cropBefore(tcut);
  cout<<"Ndata="<<data.size()<<endl;
  double t0,twidth;
  //define time range:
  double tstart,tend;
  data.getTimeLimits(tstart,tend);
  t0=data.getPeakTime();
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
      dump_view(outname,data,instate,tfinestart,tfineend,mm_samples);
    } else {
      cout<<" The -view option requires that parameters are provided."<<endl;
    }
  }

  //Set the prior:
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
  MLFitProb *llike=nullptr;
  llike = new MLFitProb(&space,&data,&prior);
  if(have_pars0){
    double ll=llike->evaluate_log(instate);
    double lp=prior.evaluate_log(instate);
    cout<<"log-Likelihood at input parameters = "<<ll<<endl;
    cout<<"log-posterior at input parameters = "<<ll+lp<<endl;
  }
  if(view)exit(0);
	  
  //Set the proposal distribution (might want to formalize some of this, now basically copied in similar form in the for the third time in a new program...)
  int Ninit;
  proposal_distribution *prop=new_proposal_distribution(Ninit,opt,&prior,&halfwidths);
  cout<<"Proposal distribution (at "<<prop<<") is:\n"<<prop->show()<<endl;

  //*************************************************
  //if(parallel_tempering)Nchain=1;

  //Prepare for chain output
  ss<<"gle_"<<outname;
  string base=ss.str();
  ss<<".dat";
  ofstream out(ss.str().c_str());
  out.precision(output_precision);

  ss.str("");ss<<base<<"_"<<Nsigma<<"_sigma.dat";
  ofstream out1sigma(ss.str().c_str());
  out1sigma.precision(output_precision);

  ss.str("");ss<<base<<"_PTstats.dat";
  ofstream outp(ss.str().c_str());
  outp.precision(output_precision);

  parallel_tempering_chains *ptc=nullptr;
  MH_chain *c=nullptr,*last_chain=nullptr;

  //*************************************************
  //Loop over Nchains
  for(int ic=0;ic<Nchain;ic++){
    proposal_distribution *chain_prop=nullptr;
    sampleable_probability_function *chain_prior=&prior;
    MLFitProb *chain_llike=nullptr;
    int chain_Nstep=Nstep,chain_Ninit=Ninit,chain_nburn=Nstep*nburn_frac;
    ofstream *chain_out=nullptr;
    
    //Run in non-progressive mode
    chain_prop=prop->clone();
    chain_prior=&prior;
    chain_llike = llike;
    chain_out=&out;	

    //Create the Chain 
    chain *cc=nullptr;
    if(parallel_tempering){
      ptc= new parallel_tempering_chains(Nptc,Tmax,swap_rate,save_every);
      if(pt_evolve_rate>0)ptc->evolve_temps(pt_evolve_rate);
      if(pt_reboot_rate>0)ptc->do_reboot(pt_reboot_rate,pt_reboot_cut,pt_reboot_thermal,pt_reboot_every,pt_reboot_grace,pt_reboot_grad,pt_reboot_blindly);
      ptc->initialize(chain_llike,chain_prior,chain_Ninit);
      //ptc->set_proposal(*chain_prop);
      cc=ptc;
      cout<<"cc=ptc="<<cc<<endl;
    } else {
      c= new MH_chain(chain_llike,chain_prior,-30,save_every);
      chain_prop->set_chain(c);
      c->initialize(chain_Ninit);
      cc=c;
    }
    cc->set_proposal(*chain_prop);
    
    cout<<"\nRunning chain "<<ic<<endl;
    //cout<<chain_llike->print_info()<<endl;
    chain_llike->reset();

    //verbose=false;
    //if(!parallel_tempering)c->reserve(chain_Nstep);//allocate
    for(int i=0;i<=chain_Nstep;i++){
      //cout<<"step="<<i<<endl;

      cc->step();
      //if(parallel_tempering)
      //ptc->step();
      //else
      //c->step(*chain_prop);
	
      //if(i>160)verbose=true;
      //cout<<"i="<<i<<"/"<<Nevery<<endl;
      if(0==i%Nevery){
	cout<<"chain "<<ic<<" step "<<i<<" MaxPosterior="<<chain_llike->best_post<<endl;
	//cout<<" capacity="<<cc->capacity()<<endl;
	//cout<<" size="<<cc->size()<<endl;
	//cout<<" i_after_burn="<<cc->i_after_burn()<<endl;
	//if(i>=30000)::verbose=true;
	/******* This isn't working right because chain-lenghts are not exactly equal to the number of steps taken.  I think that parallel-tempering swaps are recorded as 'extra' steps.  Either we decide this is wrong, or we change the use of 'i' here to something else.  Right now we assume the two are the same, resulting in significantly overlapping output series once the number of extra points exceeds Nevery **********/ 
	cc->dumpChain(*chain_out,i-Nevery+1,Nskip);
	//cc->dumpChain(cout,i-Nevery+1,Nskip);
	//if(parallel_tempering)
	//ptc->dumpChain(0,*chain_out,i-Nevery+1,Nskip);
	//else
	//c->dumpChain(*chain_out,i-Nevery+1,Nskip);
	cout<<cc->status()<<endl;
      }
    }

    if(parallel_tempering){
      //FIXME  Want to replace this with a generic chain->report() type function...
      ptc->dumpTempStats(outp);
      outp<<"\n"<<endl;
    }

    *chain_out<<"\n"<<endl;
    
    cout<<"Finished running chain "<<ic<<"."<<endl;
    //if(ic==1)::verbose=true;
    //Analysis
    //Select 1-sigma chain points
    vector<int> idx_in_Nsigma;
    ss.str("");ss<<base<<"_1_sigma_samples_"<<ic<<".dat";
    ofstream outsamples(ss.str().c_str());
    ss.str("");ss<<base<<"_1_sigma_samples_fine_"<<ic<<".dat";
    ofstream outfinesamples(ss.str().c_str());
    ss.str("");ss<<base<<"_best_"<<ic<<".dat";
    ofstream outbest(ss.str().c_str());
    ss.str("");ss<<base<<"_best_fine_"<<ic<<".dat";
    ofstream outfinebest(ss.str().c_str());

    cc->inNsigma(Nsigma,idx_in_Nsigma,chain_nburn);
    for(int i=0;i<(int)idx_in_Nsigma.size();i+=Nskip){
      int idx=idx_in_Nsigma[i];
      //cout<<"i,idx="<<i<<","<<idx<<":"<<cc->getState(idx,true).get_string()<<endl;
      //cout<<"i,idx="<<i<<","<<idx<<":"<<cc->getState(idx).get_string()<<endl;
      out1sigma<<i<<" "<<idx<<" "<<cc->getLogPost(idx,true)<<": ";
      valarray<double> p(cc->getState(idx,true).get_params());
      //valarray<double> p(cc->getState(idx).get_params());
      for(double p_j:p)out1sigma<<p_j<<" ";
      out1sigma<<endl;
    }
    out1sigma<<endl;
    outsamples.precision(output_precision);
    outfinesamples.precision(output_precision);
    double nfine=data.size()*2;
    
    for(int i=0;i<Nbest;i++){  
      int idx=idx_in_Nsigma[(rand()*idx_in_Nsigma.size())/RAND_MAX];
      state st=cc->getState(idx,true);
      //state st=cc->getState(idx);
      vector<double> p=st.get_params_vector();
      outsamples<<"#"<<st.get_string()<<endl;
      data.write(outsamples,p,true);
      outsamples<<endl;
      outfinesamples<<"#"<<st.get_string()<<endl;
      data.write(outfinesamples,p,true,nfine,tfinestart,tfineend);
      outfinesamples<<endl;
    }

    outbest.precision(output_precision);    
    outfinebest.precision(output_precision);    
    int idx=idx_in_Nsigma[0];
    state st=cc->getState(idx,true);
    //state st=cc->getState(idx);
    vector<double> p=st.get_params_vector();
    outbest<<"#"<<st.get_string()<<endl;
    data.write(outbest,p,true);
    outbest<<endl;
    outfinebest<<"#"<<st.get_string()<<endl;
    data.write(outfinebest,p,true,nfine,tfinestart,tfineend);
    outfinebest<<endl;
    
    cout<<"chain "<<ic<<": best_post "<<chain_llike->best_post<<", state="<<llike->best.get_string()<<endl;

    //clean up
    if(parallel_tempering)delete ptc;
    else delete c;
  }
  
  //oute.close();
  out.close();  
  
  //Dump summary info
  cout<<"best_post "<<llike->best_post<<", state="<<llike->best.get_string()<<endl;
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

///Set the proposal distribution. Calling routing responsible for deleting.
///Also returns choice of Ninit in first arg.
///This should eventually go to ptmcmc driver class...
proposal_distribution* new_proposal_distribution(int &Ninit, const Options &opt, sampleable_probability_function * prior, const valarray<double>*halfwidths){
  int proposal_option,SpecNinit;
  double tmixfac,reduce_gamma_by,de_eps;
  bool de_mixing=false;
  istringstream(opt.value("prop"))>>proposal_option;
  istringstream(opt.value("de_ni"))>>SpecNinit;
  istringstream(opt.value("de_eps"))>>de_eps;
  istringstream(opt.value("de_reduce_gamma"))>>reduce_gamma_by;
  istringstream(opt.value("de_Tmix"))>>tmixfac;
  de_mixing=opt.set("de_mixing");
  valarray<double> sigmas;
  if(halfwidths!=nullptr)sigmas=*halfwidths;
  else if(proposal_option<2){
    cout<<"new_proposal_distribution: Called without defining haflwidths. Cannot apply proposal option 0 or 1."<<endl;
    exit(1);
  }
  Ninit = 1;
  proposal_distribution* prop=nullptr;
  
  //if(parallel_tempering)reduce_gamma_by*=2.0;//Turned this off since mlchain.
  switch(proposal_option){
  case 0:  //Draw from prior distribution   
    cout<<"Selected draw-from-prior proposal option"<<endl;
    prop=new draw_from_dist(*prior);
    break;
  case 1:  //gaussian   
    cout<<"Selected Gaussian proposal option"<<endl;
    prop=new gaussian_prop(sigmas/8.);
    break;
  case 2:  {  //range of gaussians
    int Nprop_set=4;
    cout<<"Selected set of Gaussian proposals option"<<endl;
    vector<proposal_distribution*> gaussN(Nprop_set,nullptr);
    vector<double>shares(Nprop_set);
    double fac=1;
    for(int i=0;i<Nprop_set;i++){
      fac*=2;
      gaussN[i]=new gaussian_prop(sigmas/fac);
      shares[i]=fac;
      cout<<"  sigma="<<sigmas[0]/fac<<", weight="<<fac<<endl;
    }
    prop=new proposal_distribution_set(gaussN,shares);
    break;
  }
  case 3:{
    cout<<"Selected differential evolution proposal option"<<endl;
    //differential_evolution(bool snooker=false, double gamma_one_frac=0.1,double b_small=0.0001,double ignore_frac=0.3):snooker(snooker),gamma_one_frac(gamma_one_frac),b_small(b_small),ignore_frac(ignore_frac)
    //prop=new differential_evolution();
    differential_evolution *de=new differential_evolution(0.0,0.3,de_eps,0.0);
    de->reduce_gamma(reduce_gamma_by);
    if(de_mixing)de->support_mixing(true);
    de->mix_temperatures_more(tmixfac);
    prop=de;
    
    Ninit=SpecNinit*Npar;//Need a huge pop of samples to avoid getting stuck in a peak unless occasionally drawing from prior.
    break;
  }
  case 4:{
    cout<<"Selected differential evolution with snooker updates proposal option"<<endl;
    //differential_evolution(bool snooker=false, double gamma_one_frac=0.1,double b_small=0.0001,double ignore_frac=0.3):snooker(snooker),gamma_one_frac(gamma_one_frac),b_small(b_small),ignore_frac(ignore_frac)
    differential_evolution *de=new differential_evolution(0.1,0.3,de_eps,0.0);
    //differential_evolution *de=new differential_evolution(0.1,0.3,0.0);
    de->reduce_gamma(reduce_gamma_by);
    if(de_mixing)de->support_mixing(true);
    de->mix_temperatures_more(tmixfac);
    prop=de;
    if(false){
      vector<proposal_distribution*>props(1);
      vector<double>shares(1);
      props[0]=de;shares[0]=1.0;
      prop=new proposal_distribution_set(props,shares);    
    }
    Ninit=SpecNinit*Npar;
    break;
  }
  case 5:{
    cout<<"Selected differential evolution proposal with prior draws option"<<endl;
    //differential_evolution(bool snooker=false, double gamma_one_frac=0.1,double b_small=0.0001,double ignore_frac=0.3):snooker(snooker),gamma_one_frac(gamma_one_frac),b_small(b_small),ignore_frac(ignore_frac)
    //prop=new differential_evolution();
    vector<proposal_distribution*>props(2);
    vector<double>shares(2);
    props[0]=new draw_from_dist(*prior);
    shares[0]=0.1;
    differential_evolution *de=new differential_evolution(0.0,0.3,de_eps,0.0);
    de->reduce_gamma(reduce_gamma_by);
    if(de_mixing)de->support_mixing(true);
    de->mix_temperatures_more(tmixfac);
    cout<<"de="<<de<<endl;
    props[1]=de;
    shares[1]=0.9;
    prop=new proposal_distribution_set(props,shares);    
    Ninit=SpecNinit*Npar;
    break;
  }
  case 6:{
    cout<<"Selected differential evolution (with snooker updates) proposal with prior draws option"<<endl;
    //differential_evolution(bool snooker=false, double gamma_one_frac=0.1,double b_small=0.0001,double ignore_frac=0.3):snooker(snooker),gamma_one_frac(gamma_one_frac),b_small(b_small),ignore_frac(ignore_frac)
    //prop=new differential_evolution();
    vector<proposal_distribution*>props(2);
    vector<double>shares(2);
    props[0]=new draw_from_dist(*prior);
    shares[0]=0.1;
    differential_evolution *de=new differential_evolution(0.1,0.3,de_eps,0.0);
    de->reduce_gamma(reduce_gamma_by);
    if(de_mixing)de->support_mixing(true);
    de->mix_temperatures_more(tmixfac);
    props[1]=de;
    shares[1]=0.9;
    prop=new proposal_distribution_set(props,shares);    
    Ninit=SpecNinit*Npar;
    break;
  }
  default:
    cout<<"new_proposal_distribution: Unrecognized value: proposal_option="<<proposal_option<<endl;
    exit(1);
  }
  return prop;
}
