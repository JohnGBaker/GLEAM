#!/opt/local/bin/python
import sys
import os
import numpy as np
import math
import subprocess
import shlex
import argparse 
import time

#functions:
#read test info from file and do runs
def make_runs(fname,commandbase,tag,postprocess_only):
    filenames=[];
    lines=[];
    with open(fname) as f:
        for line in f:
            if(not line.startswith("#")):
                #print "line is: '"+line+"'"
                lines.append(line)
                cols=line.split()
                if(len(cols)<4):
                    continue
                center=int(cols[0])
                width=float(cols[1])
                qval=float(cols[2])
                Lval=float(cols[3])
                outname=tag+"_"+str(center)+"_"+str(width)+"_"+str(qval)+"_"+str(Lval)
                if(not postprocess_only):
                    command=commandbase+" -mm_center="+str(center)+" -mm_log_q="+str(qval)+" -mm_log_L="+str(Lval)+" -mm_width="+str(width)+" "+str(outname)
                    print "\n**********************************************"
                    print command
                    sys.stdout.flush()
                    subprocess.call(command.split())
                filenames.append(outname)
    return filenames,lines
# make_runs

#compare runs a multiple relative precision levels.
def compare_runs_rel(filesA,filesB,precRmax,lines):
    triesA=filesA
    triesB=filesB
    triesL=lines
    passlevels={line:0 for line in lines}
    while prec<1:
        result,fails,failsA,failsB=compare_runs(triesA,triesB,1,prec,triesL)
        print "  prec="+prec+("OK" if result else "FAIL")
        for line in triesL:
            if(not line in fails):
                passlevel[line]=prec;
        prec=prec*10.0;
        triesL=fails
        triesA=failsA
        triesB=failsB
    return passlevels   

#compare runs a multiple absolute precision levels.
def compare_runs_rel(filesA,filesB,precAmax,lines):
    triesA=filesA
    triesB=filesB
    triesL=lines
    passlevels={line:0 for line in lines}
    while prec<1:
        result,fails,failsA,failsB=compare_runs(triesA,triesB,prec,1,triesL)
        print "  prec="+prec+("OK" if result else "FAIL")
        for line in triesL:
            if(not line in fails):
                passlevel[line]=prec;
        prec=prec*10.0;
        triesL=fails
        triesA=failsA
        triesB=failsB
    return passlevels

#run ndiff on the relevant runs and collect results
def compare_runs_ndiff(filesA,filesB,lines):
    suf="_mmap.dat"
    good=True
    absrel=[]
    for i in range(len(filesA)):
        result=""
        dataA=np.loadtxt(filesA[i]+suf).flatten();
        dataB=np.loadtxt(filesB[i]+suf).flatten();
        deltas=np.array([[abs(a-b),abs(a-b)/(abs(a)+abs(b))] for a,b in zip(dataA,dataB)])
        #print deltas[:,0]
        #print deltas[:,1]
        absreli=[max(deltas[:,0]),max(deltas[:,1])]
        print "Compared files "+filesA[i]+suf+" and "+filesB[i]+suf+" : "+str(absreli) 
        absrel.append(absreli);
    return absrel
# compare_runs


#run numdiff on the relevant runs and collect results
def compare_runs(filesA,filesB,precA,precR,lines):
    suf="_mmap.dat"
    good=True
    failing_cases=[]
    failing_filesA=[]
    failing_filesB=[]
    #ngood=2
    for i in range(len(filesA)):
        result=""
        command="numdiff -H --strict -q -a "+str(precA)+" -r "+str(precR)+" "+filesA[i]+suf+" "+filesB[i]+suf
        print
        print command
        #result=subprocess.check_output(command.split())
        #nlines=result.count("\n")
        #print nlines
        #print result
        #bad=nlines>good
        bad=not 0==subprocess.call(command.split())
        #bad=False
        if (bad):
            print "FAILED"
            good=False
            failing_cases.append(lines[i])
            failing_filesA.append(filesA[i])
            failing_filesB.append(filesB[i])
        else :
            print "PASSED"
        sys.stdout.flush()
    return good, failing_cases, failing_filesA, failing_filesB
# compare_runs

##########
# setup stuff
# determine basename from first argument.
scriptpath=os.path.dirname(os.path.realpath(__file__))
parser = argparse.ArgumentParser(description='Run tests on gleam lens inversion.')
parser.add_argument('testlist',
                   help='filename for the test list')
parser.add_argument('prec', type=float,
                   help='log precision level of test')
parser.add_argument('--tag', default="",
                   help='an extra tag for the output')
parser.add_argument('--post', action='store_true',
                   help='run only post-processing')
parser.add_argument('--exec', default=scriptpath+"/../gleam",dest='executable',
                    help='path to the executable (default to expected location relative script)')
parser.add_argument('--exec-arg', default="",dest='exarg',
                    help='additional arguments for the gleam executable')
parser.add_argument('--no-quad', action='store_false',dest='do_quad',
                    help='don\'t do quad tests')

args = parser.parse_args()

#if(len(sys.argv)<3 or len(sys.argv)>4):
#    print "Usage:\n"+sys.argv[0]+" testlistfile prec_digits [execname]"
#    exit()
#testfname = sys.argv[1]
testfname=args.testlist
outdir =testfname+args.tag+"-tests/"
reportfile =testfname+args.tag+"-report.dat"
print "testfname=",testfname
print "outdir=",outdir

#prec=10**-float(sys.argv[2])
prec=10**-float(args.prec)

postprocess_only=args.post

#get execname
#if(len(sys.argv)>3):
#    execname=sys.argv[3]
#else:
#    execname="../src/gleam"
execname=args.executable
print "execname="+execname

#create/empty output directory
if(execname.count("echo")<1 and not postprocess_only):
    command= "mkdir "+outdir
    print command
    subprocess.call(command.split())
    #os.chdir(outdir)
    for fileName in os.listdir(outdir):
        os.remove(outdir+fileName)


#process different types of runs
command= execname+" "+args.exarg+" -magmap -poly=true -precision=16 -GLB_rWide=5 "
tag="poly_r5.0"
start = time.time()
poly5files,lines=make_runs(testfname,command,outdir+tag,postprocess_only)
end=time.time()
timeP5=end-start;

command= execname+" "+args.exarg+" -magmap -poly=true -precision=16 -GLB_rWide=4.5 "
tag="poly_r4.5"
start = time.time()
poly4files,lines=make_runs(testfname,command,outdir+tag,postprocess_only)
timeP4=time.time()-start;

if(args.do_quad):
    command= execname+"_quad "+args.exarg+" -magmap -poly=true -precision=16 -GLB_rWide=5 "
    tag="qpoly_r5.0"
    start = time.time()
    qpoly5files,lines=make_runs(testfname,command,outdir+tag,postprocess_only)
    timeQP5=time.time()-start;

#command= execname+"_quad "+args.exarg+" -magmap -poly=true -precision=16 -GLB_rWide=4.5 "
#tag="qpoly_r4.5"
#start = time.time()
#qpoly4files,lines=make_runs(testfname,command,outdir+tag,postprocess_only)
#timeQP4=time.time()-start;

command= execname+" "+args.exarg+" -magmap -precision=16 -GLB_rWide=5 "
tag="int_r5.0"
start = time.time()
int5files,lines=make_runs(testfname,command,outdir+tag,postprocess_only)
timeI5=time.time()-start;

#compare
good=goodr=True

with open(reportfile,'w') as report:
    print" Comparing WideBinary versus WittMao."
    report.write("\n# Comparing WideBinary versus WittMao.\n")
    resultWBWM=compare_runs_ndiff(poly5files,poly4files,lines)
    for l,r in zip(lines,resultWBWM):
        ls=l.split()
        report.write(ls[0]+" "+ls[1]+" "+ls[2]+" "+ls[3]+" "+str(r[0])+" "+str(r[1])+"\n")
    print"\n Comparing Polynomial versus Integration"
    report.write("\n# Comparing Polynomial versus Integration")
    resultPI=compare_runs_ndiff(poly5files,int5files,lines)
    for l,r in zip(lines,resultPI):
        ls=l.split()
        report.write(ls[0]+" "+ls[1]+" "+ls[2]+" "+ls[3]+" "+str(r[0])+" "+str(r[1])+"\n")
    if(args.do_quad):
        print" Comparing quad versus double precision."
        report.write("\n# Comparing quad versus double precision.")
        resultQD=compare_runs_ndiff(poly5files,qpoly5files,lines)
        for l,r in zip(lines,resultQD):
            ls=l.split()
            report.write(ls[0]+" "+ls[1]+" "+ls[2]+" "+ls[3]+" "+str(r[0])+" "+str(r[1])+"\n")
        print" Comparing Quad-precision Polynomial versus Integration"
        report.write("\n# Comparing Quad-precision Polynomial versus Integration")
        resultQPI = compare_runs_ndiff(qpoly5files,int5files,lines)
        for l,r in zip(lines,resultQPI):
            ls=l.split()
            report.write(ls[0]+" "+ls[1]+" "+ls[2]+" "+ls[3]+" "+str(r[0])+" "+str(r[1])+"\n")

if( not postprocess_only):
    print "\nTiming summary:"
    print "         Poly 5.0 - ",timeP5
    print "         Poly 4.5 - ",timeP4
    if(args.do_quad):
        print "    Quad Poly 5.0 - ",timeQP5
#    print "    Quad Poly 4.5 - ",timeQP4
    print "          Int 5.0 - ",timeI5


print "\nTest summary at relative prec="+str(prec)+":"
print "         WB vs WM: "+("OK "+str(max(np.array(resultWBWM)[:,1])) if max(np.array(resultWBWM)[:,1])<=prec else "FAIL")
for l,r in zip(lines,resultWBWM):
    ls=l.split()
    if(r[1]>prec):
        print ls[0]+" "+ls[1]+" "+ls[2]+" "+ls[3]+" "+str(r[1])
print "      poly vs int: "+("OK "+str(max(np.array(resultPI)[:,1])) if max(np.array(resultPI)[:,1])<=prec else "FAIL")
for l,r in zip(lines,resultPI):
    ls=l.split()
    if(r[1]>prec):
        print ls[0]+" "+ls[1]+" "+ls[2]+" "+ls[3]+" "+str(r[1])
if(args.do_quad):
    print "   quad vs double: "+("OK "+str(max(np.array(resultQD)[:,1])) if max(np.array(resultQD)[:,1])<=prec else "FAIL")
    for l,r in zip(lines,resultQD):
        ls=l.split()
        if(r[1]>prec):
            print ls[0]+" "+ls[1]+" "+ls[2]+" "+ls[3]+" "+str(r[1])
    print " quad poly vs int: "+("OK "+str(max(np.array(resultQPI)[:,1])) if  max(np.array(resultQPI)[:,1])<=prec else "FAIL")
    for l,r in zip(lines,resultQPI):
        ls=l.split()
        if(r[1]>prec):
            print ls[0]+" "+ls[1]+" "+ls[2]+" "+ls[3]+" "+str(r[1])


print "\nTest summary at absolute prec="+str(prec)+":"
print "         WB vs WM: "+("OK "+str(max(np.array(resultWBWM)[:,0])) if max(np.array(resultWBWM)[:,0])<=prec else "FAIL")
for l,r in zip(lines,resultWBWM):
    ls=l.split()
    if(r[0]>prec):
        print ls[0]+" "+ls[1]+" "+ls[2]+" "+ls[3]+" "+str(r[0])
print "      poly vs int: "+("OK "+str(max(np.array(resultPI)[:,0])) if max(np.array(resultPI)[:,0])<=prec else "FAIL")
for l,r in zip(lines,resultPI):
    ls=l.split()
    if(r[0]>prec):
        print ls[0]+" "+ls[1]+" "+ls[2]+" "+ls[3]+" "+str(r[0])
if(args.do_quad):
    print "   quad vs double: "+("OK "+str(max(np.array(resultQD)[:,0])) if max(np.array(resultQD)[:,0])<=prec else "FAIL")
    for l,r in zip(lines,resultQD):
        ls=l.split()
        if(r[0]>prec):
            print ls[0]+" "+ls[1]+" "+ls[2]+" "+ls[3]+" "+str(r[0])
    print " quad poly vs int: "+("OK "+str(max(np.array(resultQPI)[:,0])) if  max(np.array(resultQPI)[:,0])<=prec else "FAIL")
    for l,r in zip(lines,resultQPI):
        ls=l.split()
        if(r[0]>prec):
            print ls[0]+" "+ls[1]+" "+ls[2]+" "+ls[3]+" "+str(r[0])



