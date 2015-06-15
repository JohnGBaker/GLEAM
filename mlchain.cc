
#include "mlfit.hh"
#include <valarray>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>

#ifdef USE_PTMCMCS
#include "ptMCMCs.hh"
#else
#include "mcmc.hh"
#include "proposal_distribution.hh"
#include "chain.hh"
#endif


//#include <mcheck.h>

using namespace std;

typedef initializer_list<double> dlist;
typedef initializer_list<int> ilist;

bool debug = false;
bool debugint = false;

shared_ptr<Random> globalRNG;//used for some debugging... 

//Control parameters (see more in main() near line 110)
const int Npar=9;
double nburn_frac=0.500;
int Nstep; //Number of steps attemped in each chain
int Nprop_set=4;
double narrow_by=2;
double reduce_gamma_by=0.5;
//double reduce_gamma_by=1.0;
//const int SpecNinit=100;
//const int SpecNinit=10;
const int SpecNinit=10;
const bool use_remapped_r0=true;
const bool use_remapped_q=true;
const double q0=1e4;

int Nlead_args=5;

bool verbose=false;//probably not needed.

  //Some initializing options
#ifdef PARALLEL_TEMPERING    
const int Nptc=16; //Number of parallel temperature chains
const int Tmax=100; //maximum temperature
//const int Nptc=64; //Number of parallel temperature chains
//const int Tmax=400; //maximum temperature
//const int Nptc=64; //Number of parallel temperature chains
//const int Tmax=5000; //maximum temperature
const int parallel_tempering=true;
int Nchain=3;
#else
const int Nptc=1;
const int Tmax=0;
const bool parallel_tempering=false;
int Nchain=6;
#endif

int Nevery=500;
int idump_sampling=10;
int Nsigma=1;
int Nbest=10;


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
    //cout<<"Exaluating likelihood (this="<<this<<") for params"<<s.get_string()<<endl;
    valarray<double>params=s.get_params();
    bool integrate=false;
#ifdef INTEGRATE
    integrate=true;
#endif
    clock_t tstart=clock();
    double result=data->evaluate_log(params,integrate);
    double post=result;
    if(prior)post+=prior->evaluate_log(s);
    clock_t tend=clock();
    double eval_time = (tend-tstart)/(double)CLOCKS_PER_SEC;
    //#pragma omp critical
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
  double seed;
  string datafile;
  int proposal_option;
  const int NparRead=Npar; 

  int expectedNargs=1+Nlead_args+NparRead;
  if((argc!=expectedNargs&&argc!=expectedNargs-NparRead)||!(istringstream(argv[2])>>seed)){
    cout<<"argc="<<argc<<" (expected "<<expectedNargs<<" or "<<expectedNargs-NparRead<<")"<<endl;
    if(argc>2)
      cout<<"argv[1]="<<argv[1]<<endl;
    cout<<"Usage: \n"<<argv[0]<<" datafile seed N-steps proposal_option outname [ I0 Fs Fn logq logL r0 phi tE tpass ]";
    cout<<"\nseed should bin in [0,1)."<<endl;
    return -1;
  }
  //lead args
  datafile=argv[1];
  istringstream(argv[3])>>Nstep;
  //int nburn=Nstep*nburn_frac;
  istringstream(argv[4])>>proposal_option;
  string outname=argv[5];
  //params
  valarray<double>params(Npar);
  //set pointers to param values (might not need)
  //double * par_I0=&params[0],*par_Fs=&params[1],*par_Fn=&params[2],*par_logq=&params[3],*par_logL=&params[4],*par_r0=&params[5],*par_phi=&params[6],*par_tE=&params[7],*par_tpass=&params[8];

  bool have_pars0=false;
  if(argc>1+Nlead_args){
    for(int i=0;i<NparRead;i++)stringstream(argv[i+1+Nlead_args])>>params[i];
    have_pars0=true;
  }
    
  cout.precision(15);
  cout<<"\noutname = '"<<outname<<"'"<<endl;
  cout<<"seed="<<seed<<endl;
  //cout<<"I0  = "<<*par_I0<<"\nFs  = "<<*par_Fs<<"\nlogq   = "<<*par_logq<<"\nlogL   = "<<*par_logL<<"\nr0  = "<<*par_r0<<"\nphi = "<<*par_phi<<"\ntE  = "<<*par_tE<<"\ntpas= "<<*par_tpass<<endl;
  
  
  ProbabilityDist::setSeed(seed);
  globalRNG.reset(ProbabilityDist::getPRNG());//just for safety to keep us from deleting main RNG in debugging.
  cout<<"globalRNG="<<globalRNG.get()<<endl;
  
  //Load the data
  OGLEdata data(datafile);
  cout<<"&data="<<&data<<endl;
  double t0,twidth,r0s=6.0;
  //define time range:
  double tstart,tend;
  data.getTimeLimits(tstart,tend);
  t0=data.getPeakTime();
  double tfinestart=t0-tstart+tend;
  double tfineend=t0+tstart-tend; //tfine range is twice data range centered on t0
  twidth=10;//look mostly within 10 days of peak;
  if(use_remapped_r0){
    data.remap_r0(2.0);
    r0s=1.0;
  }
  if(use_remapped_q){
    data.remap_q(q0);
  }
  //Set up the parameter space
  stateSpace space(Npar);
  {
    string names[]={"I0","Fs","Fn","logq","logL","r0","phi","tE","tpass"};
    if(use_remapped_r0)names[5]="s(r0)";
    if(use_remapped_q)names[3]="s(1+q)";    
    space.set_names(names);  
    space.set_bound(6,boundary(boundary::wrap,boundary::wrap,0,2*M_PI));//set 2-pi-wrapped space for phi.
  }
  cout<<"&space="<<&space<<endl; 
  cout<<"Parameter space:\n"<<space.show()<<endl;

  //Initial parameter state
  state instate(&space,params);
  if(have_pars0){
    cout<<"Initial parameters:"<<instate.show()<<endl;
    cout<<" ... are basically ignored..."<<endl;
  }

  //Set the prior:
  const int uni=mixed_dist_product::uniform, gauss=mixed_dist_product::gaussian, pol=mixed_dist_product::polar;
  //                                     I0      Fs     Fn  logq*   logL      r0     phi     tE   tpass  
  valarray<double>    centers((dlist){ 18.0,    0.5,   0.5,   0.0,   0.0,r0s/2.0,   M_PI,  50.0,  t0     });
  valarray<double> halfwidths((dlist){  5.0,    0.5,   0.5,   3.0,   2.0,r0s/2.0,   M_PI,  50.0,  twidth });
  valarray<int>         types((ilist){gauss,    uni,   uni, gauss, gauss,    uni,    uni,   uni,  gauss  });
  if(use_remapped_q){//* 
    double qq=2.0/(q0+1.0);
    double ds=0.5/(1.0+qq*qq); //ds=(1-s(q=1))/2
    centers[3]=1.0-ds;
    halfwidths[3]=ds;          //ie range=[s(q=1),s(q=inf)=1.0]
    types[3]=uni;
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

  //Set the proposal distribution (might want to formalize some of this, now basically copied in similar form in the for the third time in a new program...)
  valarray<double> sigmas(halfwidths);
  int Ninit = 1;
  proposal_distribution* prop=nullptr;
  if(parallel_tempering)reduce_gamma_by*=2.0;//Don't know why but acceptance ratio come out smaller with pt
  switch(proposal_option){
  case 0:  //Draw from prior distribution   
    cout<<"Selected draw-from-prior proposal option"<<endl;
    prop=new draw_from_dist(prior);
    break;
  case 1:  //gaussian   
    cout<<"Selected Gaussian proposal option"<<endl;
    prop=new gaussian_prop(sigmas/8.);
    break;
  case 2:  {  //range of gaussians
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
    differential_evolution *de=new differential_evolution(0.0,0.3,0.0001,0.0);
    de->reduce_gamma(reduce_gamma_by);
    prop=de;
    
    Ninit=SpecNinit*Npar;//Need a huge pop of samples to avoid getting stuck in a peak unless occasionally drawing from prior.
    break;
  }
  case 4:{
    cout<<"Selected differential evolution with snooker updates proposal option"<<endl;
    //differential_evolution(bool snooker=false, double gamma_one_frac=0.1,double b_small=0.0001,double ignore_frac=0.3):snooker(snooker),gamma_one_frac(gamma_one_frac),b_small(b_small),ignore_frac(ignore_frac)
    differential_evolution *de=new differential_evolution(0.1,0.3,0.0001,0.0);
    //differential_evolution *de=new differential_evolution(0.1,0.3,0.0);
    de->reduce_gamma(reduce_gamma_by);
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
    props[0]=new draw_from_dist(prior);
    shares[0]=0.1;
    differential_evolution *de=new differential_evolution(0.0,0.3,0.0001,0.0);
    de->reduce_gamma(reduce_gamma_by);
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
    props[0]=new draw_from_dist(prior);
    shares[0]=0.1;
    differential_evolution *de=new differential_evolution(0.1,0.3,0.0001,0.0);
    de->reduce_gamma(reduce_gamma_by);
    props[1]=de;
    shares[1]=0.9;
    prop=new proposal_distribution_set(props,shares);    
    Ninit=SpecNinit*Npar;
    break;
  }
  default:
    cout<<"Unrecognized value: proposal_option="<<proposal_option<<endl;
    exit(1);
  }
  cout<<"Proposal distribution (at "<<prop<<") is:\n"<<prop->show()<<endl;

  //*************************************************
  //if(parallel_tempering)Nchain=1;

  //Prepare for chain output
  ostringstream ss("");
  if(parallel_tempering)ss<<"pt";
  ss<<"mlc_"<<outname;
  string base=ss.str();
  ss<<".dat";
  ofstream out(ss.str().c_str());
  out.precision(13);

  ss.str("");ss<<base<<"_"<<Nsigma<<"_sigma.dat";
  ofstream out1sigma(ss.str().c_str());
  out1sigma.precision(13);

#ifdef PARALLEL_TEMPERING    
  ss.str("");ss<<base<<"_PTstats.dat";
  ofstream outp(ss.str().c_str());
  outp.precision(13);
#endif

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
      ptc= new parallel_tempering_chains(Nptc,Tmax);
      ptc->initialize(chain_llike,chain_prior,chain_Ninit);
      //ptc->set_proposal(*chain_prop);
      cc=ptc;
      cout<<"cc=ptc="<<cc<<endl;
    } else {
      c= new MH_chain(chain_llike,chain_prior);
      chain_prop->set_chain(c);
      c->initialize(chain_Ninit);
      cc=c;
    }
    cc->set_proposal(*chain_prop);
    
    cout<<"\nRunning chain "<<ic<<endl;
    //cout<<chain_llike->print_info()<<endl;
    chain_llike->reset();

    verbose=false;
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
	cc->dumpChain(*chain_out,i-Nevery+1,idump_sampling);
	//cc->dumpChain(cout,i-Nevery+1,idump_sampling);
	//if(parallel_tempering)
	//ptc->dumpChain(0,*chain_out,i-Nevery+1,idump_sampling);
	//else
	//c->dumpChain(*chain_out,i-Nevery+1,idump_sampling);
	cout<<cc->status()<<endl;
      }
    }
#ifdef PARALLEL_TEMPERING
    if(parallel_tempering){
      ptc->dumpTempStats(outp);
      outp<<"\n"<<endl;
    }
#endif
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
    for(int i=0;i<(int)idx_in_Nsigma.size();i++){
      int idx=idx_in_Nsigma[i];
      out1sigma<<i<<" "<<idx<<" "<<cc->getLogPost(idx)<<": ";
      valarray<double> p(cc->getState(idx).get_params());
      for(double p_j:p)out1sigma<<p_j<<" ";
      out1sigma<<endl;
    }
    out1sigma<<endl;
    outsamples.precision(13);
    outfinesamples.precision(13);
    double nfine=data.size()*2;
    
    for(int i=0;i<Nbest;i++){  
      int idx=idx_in_Nsigma[(rand()*idx_in_Nsigma.size())/RAND_MAX];
      state st=cc->getState(idx);
      vector<double> p=st.get_params_vector();
      outsamples<<"#"<<st.get_string()<<endl;
      data.write(outsamples,p,true);
      outsamples<<endl;
      outfinesamples<<"#"<<st.get_string()<<endl;
      data.write(outfinesamples,p,true,nfine,tfinestart,tfineend);
      outfinesamples<<endl;
    }

    outbest.precision(13);    
    outfinebest.precision(13);    
    int idx=idx_in_Nsigma[0];
    state st=cc->getState(idx);
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


