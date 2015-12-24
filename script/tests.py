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
                    command=commandbase+" -mm_center="+str(center)+" "+str(outname)+" "+str(qval)+" "+str(Lval)+" "+str(width)
                    print "\n**********************************************"
                    print command
                    sys.stdout.flush()
                    subprocess.call(command.split())
                filenames.append(outname)
    return filenames,lines
# make_runs

#run numdiff on the relevant runs and collect results
def compare_runs(filesA,filesB,precA,precR,lines):
    suf="_mmap.dat"
    good=True
    failing_cases=[]
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
        else :
            print "PASSED"
        sys.stdout.flush()
    return good, failing_cases            
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

args = parser.parse_args()

#if(len(sys.argv)<3 or len(sys.argv)>4):
#    print "Usage:\n"+sys.argv[0]+" testlistfile prec_digits [execname]"
#    exit()
#testfname = sys.argv[1]
testfname=args.testlist
outdir =testfname+args.tag+"-tests/"
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
command= execname+" "+args.exarg+" -magmap -poly=true -precision=16 -mm_lens_rWB=5 "
tag="poly_r5.0"
start = time.time()
poly5files,lines=make_runs(testfname,command,outdir+tag,postprocess_only)
end=time.time()
timeP5=end-start;

command= execname+" "+args.exarg+" -magmap -poly=true -precision=16 -mm_lens_rWB=4.5 "
tag="poly_r4.5"
start = time.time()
poly4files,lines=make_runs(testfname,command,outdir+tag,postprocess_only)
timeP4=time.time()-start;

command= execname+"_quad "+args.exarg+" -magmap -poly=true -precision=16 -mm_lens_rWB=5 "
tag="qpoly_r5.0"
start = time.time()
qpoly5files,lines=make_runs(testfname,command,outdir+tag,postprocess_only)
timeQP5=time.time()-start;

#command= execname+"_quad "+args.exarg+" -magmap -poly=true -precision=16 -mm_lens_rWB=4.5 "
#tag="qpoly_r4.5"
#start = time.time()
#qpoly4files,lines=make_runs(testfname,command,outdir+tag,postprocess_only)
#timeQP4=time.time()-start;

command= execname+" "+args.exarg+" -magmap -precision=16 -mm_lens_rWB=5 "
tag="int_r5.0"
start = time.time()
int5files,lines=make_runs(testfname,command,outdir+tag,postprocess_only)
timeI5=time.time()-start;

#compare
good=goodr=True

print" Comparing WideBinary versus WittMao."
resultWBWM,failsWBWM=compare_runs(poly5files,poly4files,prec,1,lines)
#print "resultWBWM=",resultWBWM
print "  "+("OK" if resultWBWM else "FAIL")
good=good and resultWBWM
resultWBWMr,failsWBWMr=compare_runs(poly5files,poly4files,1,prec,lines)
print "  "+("OK" if resultWBWMr else "FAIL")
goodr=goodr and resultWBWMr

print" Comparing quad versus double precision."
resultQD,failsQD=compare_runs(poly5files,qpoly5files,prec,1,lines)
print "  "+("OK" if resultQD else "FAIL")
good=good and resultQD
resultQDr,failsQDr=compare_runs(poly5files,qpoly5files,1,prec,lines)
print "  "+("OK" if resultQDr else "FAIL")
goodr=goodr and resultQDr

#print" Comparing WideBinary versus WittMao at quad precision."
#resultWBWMQ,failsWBWMQ=compare_runs(qpoly5files,qpoly4files,prec,1,lines)
#print "  "+("OK" if resultWBWMQ else "FAIL")
#good=good and resultWBWMQ
#resultWBWMQr,failsWBWMQr=compare_runs(qpoly5files,qpoly4files,1,prec,lines)
#print "  "+("OK" if resultWBWMQr else "FAIL")
#goodr=goodr and resultWBWMQr

print" Comparing Polynomial versus Integration"
resultPI,failsPI=compare_runs(poly5files,int5files,prec,1,lines)
print "  "+("OK" if resultPI else "FAIL")
good=good and resultPI
resultPIr,failsPIr=compare_runs(poly5files,int5files,1,prec,lines)
print "  "+("OK" if resultPIr else "FAIL")
goodr=goodr and resultPIr

print" Comparing Quad-precision Polynomial versus Integration"
resultQPI,failsQPI=compare_runs(qpoly5files,int5files,prec,1,lines)
print "  "+("OK" if resultPI else "FAIL")
good=good and resultQPI
resultQPIr,failsQPIr=compare_runs(qpoly5files,int5files,1,prec,lines)
print "  "+("OK" if resultPIr else "FAIL")
goodr=goodr and resultQPIr

if( not postprocess_only):
    print "\nTiming summary:"
    print "         Poly 5.0 - ",timeP5
    print "         Poly 4.5 - ",timeP4
    print "    Quad Poly 5.0 - ",timeQP5
#    print "    Quad Poly 4.5 - ",timeQP4
    print "          Int 5.0 - ",timeI5


print "\nTest summary at absolute prec="+str(prec)+":"
print "         WB vs WM: "+("OK" if resultWBWM else "FAIL")
for s in failsWBWM:
    sys.stdout.write("           "+s)
print "   quad vs double: "+("OK" if resultQD else "FAIL")
for s in failsQD:
    sys.stdout.write("           "+s)
#print "  WB vs WM (quad): "+("OK" if resultWBWMQ else "FAIL")
#for s in failsWBWMQ:
#    sys.stdout.write("           "+s)
print "      poly vs int: "+("OK" if resultPI else "FAIL")
for s in failsPI:
    sys.stdout.write("           "+s)
print " quad poly vs int: "+("OK" if resultQPI else "FAIL")
for s in failsQPI:
    sys.stdout.write("           "+s)
print "Overall absolute prec test: "+("PASSED" if good else "FAILED")


print "\nTest summary at relative prec="+str(prec)+":"
print "         WB vs WM: "+("OK" if resultWBWMr else "FAIL")
for s in failsWBWMr:
    sys.stdout.write("           "+s)
print "   quad vs double: "+("OK" if resultQDr else "FAIL")
for s in failsQDr:
    sys.stdout.write("           "+s)
#print "  WB vs WM (quad): "+("OK" if resultWBWMQr else "FAIL")
#for s in failsWBWMQr:
#    sys.stdout.write("           "+s)
print "      poly vs int: "+("OK" if resultPIr else "FAIL")
for s in failsPIr:
    sys.stdout.write("           "+s)
print " quad poly vs int: "+("OK" if resultQPIr else "FAIL")
for s in failsQPIr:
    sys.stdout.write("           "+s)
print "Overall relative prec test: "+("PASSED" if goodr else "FAILED")


