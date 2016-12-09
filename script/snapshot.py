#!/opt/local/bin/python
import sys
import os
import numpy as np
import math
import subprocess
import argparse

#functions:

#read step params from 
def get_step_pars(fname,stop_size=-1):
    fsize = os.stat(fname).st_size
    print "fsize=",fsize
    if(stop_size>=0):#sample from some specifed point in the middle of the file, instead of the end
        fsize=stop_size
    #print "fsize=",fsize
    npar  = 11
    bufsize = 15*(npar+5)
    with open(fname,'rb') as f:
        #print "bufsize=",bufsize
        if bufsize > fsize:
            bufsize = fsize-1
        #print "bufsize=",bufsize
        #print "seeking=",fsize-bufsize
        data = []
        f.seek(fsize-bufsize)
        data.extend(f.readlines(bufsize))
        line=data[len(data)-3]
        print "line is:",line 
        step= int(line.split()[0])
        pars= np.array(line.split()[5:])
        #pars=pars[:-1]#last column is temp
    #print "step="+step
    return step,pars
# get_step_pars

#count the number of completed chains in the output file, based on number blank lines
def count_chains(fname):
    nblank=0
    chain_ends=[0,0,0]
    pos=0
    with open(fname,'rb') as f:
        for line in f:
            if(line.startswith("\n")):
                #print "line:",line
                chain_ends[int(nblank/2)]=pos
                #print "set: chain_ends["+str(int(nblank/2))+"]=",pos
                nblank+=1
            pos+=len(line)-1 #the "-1" is a hack because seek/tell/counting bytes all don't work together in an effective way. This does not seem to be related to the buffering issue on the tell() side with file iteration.
        #bufsize = 15*(9+5)
        #print "testseeking=",chain_ends[0]-bufsize
        #data = []
        #f.seek(chain_ends[0]-bufsize)
        #data.extend(f.readlines(bufsize))        
        #print "testline is:",data[len(data)-3] 
    print "chain_ends=",chain_ends
    return nblank/2,chain_ends; #two blanks after completed chain.
#count_chains

#read arguments from outfile 
def get_flags(fname):
    d={}
    with open(fname) as f:
        line=f.readline()
        if(not line.count("flags=")):
           print "I can't find the flags in the .out file: '",fname,"'"
           sys.exit()
        lines=f.readlines(2000)
        for line in lines:
            entries=line.split();
            if(len(entries)==0):
                break
            entry=line.split()[0] #eliminate whitespace
            vals=entry.split(":")
            if(not vals[1].count("not")>0):
                d[vals[0]]=vals[1]
                print vals[0],"->",vals[1]
    return d
# get_flags
  
##########
# setup stuff
# determine basename from first argument.

parser = argparse.ArgumentParser(description='Provide snapshot of chain state.')
parser.add_argument('fname', metavar='chain_file', type=str, 
                   help='chain file path')
parser.add_argument('datafile', metavar='data_file', type=str, nargs='?', default='',
                   help='lightcurve data file path')
parser.add_argument('-i','--ichain', metavar='chain_index', type=int, nargs=1, default=-1,
                   help='index of chain to view')
parser.add_argument('-e','--execname', metavar='executable_file', type=str, nargs=1, default="../../src/gleam/gleam",
                   help='executable file path')
parser.add_argument('-p','--points', metavar='n', type=int, default=300,
                   help='specify number of sample points')

args = parser.parse_args()
print(args)

#fname = sys.argv[1]
fname=args.fname
if(fname.endswith("_t0.out")):
    #we assume ".out does not otherwise appear in the name... that
    #could check... if(fname.count(".out")>1):...
    basename=fname.replace("_t0.out","")
    #make a guess of the chainfile name
    fname=fname.replace("_t0.out",".dat")
elif(fname.endswith(".dat")):
    #basename=fname.replace("gle_","")
    basename=fname
    basename=basename.replace("_t0.dat","")
print "basename="+basename

#get execname
#if(len(sys.argv)>2):
#    execname=sys.argv[2]
#else:
#    execname="../../src/gleam"
execname=args.execname[0]
print "execname=",execname

#get datafile
#if(len(sys.argv)>3):
#    datafile=sys.argv[3]
#else:
#    print "Need a method for automatically identifying the datafile, or provide it as third argument."
#    sys.exit()
datafile=args.datafile
print "datafile=",datafile

#get optional specification of which chain to take
#if(len(sys.argv)>4):
#    ichain=int(sys.argv[4])
#else:
#    ichain=-1
ichain=args.ichain[0]
print "ichain="+str(ichain)

npoints=args.points
print "npoints=",npoints

#get gnuplot script path
gpscript=sys.argv[0].replace(".py",".gp")

#determine arguments
d=get_flags(basename+".out")
eargs=""
for key in d.keys():
    eargs+=" -"+key+"="+d[key]
print "eargs=",eargs

#get pars and other info from .dat file
cc,cends=count_chains(fname)
print "cends=",cends," cc=",cc
iend=-1
if(ichain>=0 and cc>ichain):
    iend=cends[ichain]
    cc=ichain
print "iend=",iend," cc=",cc

step,pars = get_step_pars(fname,iend)
resultname=basename+"_c"+str(cc)+"_"+str(int(math.floor(step/1000)))+"k"
parstr=""
for val in pars:
    parstr+=" "+str(val)
parfile=resultname+".pars"
with open(parfile,"w") as f:f.write(parstr+"\n")

print "step 0"

command= execname+" -view -mm_samples="+str(npoints)+" -stateFile="+parfile+" "+eargs+" "+datafile+" "+resultname
print command
print command.split()
sys.stdout.flush()
subprocess.call(command.split())
print "step 1"


command="gnuplot -e basename='"+resultname+"' "+gpscript
print command
print command.split()
sys.stdout.flush()
subprocess.call(command.split())
print "step 2"

command="open "+resultname+".pdf"
print command
print command.split()
sys.stdout.flush()
subprocess.call(command.split())
