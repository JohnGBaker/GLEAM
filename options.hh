//Code for parsing command-line options.
//Written by John Baker/NASA-GSFC (2010-2014)
#include <string>
#include <iostream>
#include <sstream>
#include <map>
using namespace std;


class Options;
class Option {
  friend class Options;
 private:
  string name;
  string info;
  string value;
  bool have_default;
  bool is_set;
 public:
  Option();
  Option(const char *name, const char *info, const char* vdefault="<no default>");
};

Option::Option():name(""),info(""),value(""),have_default(false),is_set(false){};

Option::Option(const char *name, const char *info, const char *vdefault):name(name),info(info),value(vdefault){
  is_set= false;
  have_default=false;
  if(string(vdefault).compare("<no default>")!=0)have_default=true;
}

class Options {
 private:
  map<string,Option> flags;
 public:
  Options(){};
  void add(Option opt);
  bool parse(int & argc, char* argv[]);
  bool set(const string &name, string &return_value)const;
  bool set(const string &name)const{
    string dummy;
    return set(name,dummy);
  };
  string value(const string &name)const{
    string s("");
    set(name,s);
    return s;
  };
  string print_usage()const;
  string report();
};

void Options::add(Option opt){
  flags[ opt.name ] = opt;
}

bool Options::set(const string &name, string & value)const{
  const Option *opt =&flags.find(name)->second;
  if(flags.count(name)==0){
    cerr<<"Options: Error no option '"<<name<<"'."<<endl;
    return false;
  }
  //opt=
  if(opt->is_set||opt->have_default){
    value=opt->value;
    return true;
  }
  return false;
}

string Options::report(){
  ostringstream s;
  for (auto it:flags)
    s<< " " << it.first << ':' << (it.second.is_set?it.second.value:"(not set)") << '\n';
  return s.str();
}

bool Options::parse(int & argc, char* argv[]){
  bool fail=false;
  int count=0;
  for(int i=1;i<argc;i++){
    if(argv[i][0]!='-')break;
    count++;
    string flag( & argv[i][1] );
    unsigned int pos=flag.find_first_of("=",1);
    string name=flag.substr(0,pos);
    if(flags.count(name)==0){
      cerr<<"Option '"<<name<<"' not recognized."<<endl;
      fail=true;
    } else {
      Option *opt=&flags[name];
      opt->is_set=true;
      if(pos!=(unsigned int)string::npos){
	//cout<<"pos="<<pos<<endl;
	opt->value=flag.substr(pos+1);
	//cout<<"value='"<<opt->value<<"'"<<endl;
      }
      else
	opt->value="true";
    }
  }
  argc-=count;
  for(int i=0;i<argc;i++)argv[i+1]=argv[i+1+count];
  //cout<<"Counted "<<count<<" flags."<<endl;
  //cout<<"And "<<argc<<" arguments."<<endl;
  return fail;
}

string Options::print_usage()const{
  ostringstream os("");
  os<<"Options:\n";
  map<string,Option>::const_iterator i;

  for(i=flags.begin();i!=flags.end();i++){
    os<<"  -"<<(*i).second.name<<":\n";
    os<<"\t"<<(*i).second.info<<":\n";
  }
  return os.str();
}

///This is an interface for classes (objects?) to provide options that can be realized (e.g.) on the command line.
class Optioned{
public:
  ///Add options to the program's option list.
  virtual void addOptions(Options &opt)=0;
  virtual void processOptions(const Options &opt)=0;
};

