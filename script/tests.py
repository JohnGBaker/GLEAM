#!/opt/local/bin/python
import sys
import os
import numpy as np
import math
import subprocess
import shlex

#functions:
#read test info from file and do runs
def make_runs(fname,commandbase,tag):
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
        command="numdiff --strict -q -a "+str(precA)+" -r "+str(precR)+" "+filesA[i]+suf+" "+filesB[i]+suf
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
if(len(sys.argv)<3 or len(sys.argv)>4):
    print "Usage:\n"+sys.argv[0]+" testlistfile prec_digits [execname]"
    exit()
testfname = sys.argv[1]
outdir =testfname+"-tests/"

print "testfname=",testfname
print "outdir=",outdir

prec=10**-float(sys.argv[2])

#get execname
if(len(sys.argv)>3):
    execname=sys.argv[3]
else:
    execname="../src/gleam"
print "execname="+execname

#create/empty output directory
if(execname.count("echo")<1):
    command= "mkdir "+outdir
    print command
    subprocess.call(command.split())
    #os.chdir(outdir)
    for fileName in os.listdir(outdir):
        os.remove(outdir+fileName)


#process different types of runs

command= execname+" -magmap -poly=true -precision=16 -mm_lens_rWB=5 "
tag="poly_r5.0"
poly5files,lines=make_runs(testfname,command,outdir+tag)

command= execname+" -magmap -poly=true -precision=16 -mm_lens_rWB=4.5 "
tag="poly_r4.5"
poly4files,lines=make_runs(testfname,command,outdir+tag)

command= execname+"_quad -magmap -poly=true -precision=16 -mm_lens_rWB=5 "
tag="qpoly_r5.0"
qpoly5files,lines=make_runs(testfname,command,outdir+tag)

command= execname+"_quad -magmap -poly=true -precision=16 -mm_lens_rWB=4.5 "
tag="qpoly_r4.5"
qpoly4files,lines=make_runs(testfname,command,outdir+tag)

command= execname+" -magmap -precision=16 -mm_lens_rWB=5 "
tag="int_r5.0"
int5files,lines=make_runs(testfname,command,outdir+tag)

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

print" Comparing WideBinary versus WittMao at quad precision."
resultWBWMQ,failsWBWMQ=compare_runs(qpoly5files,qpoly4files,prec,1,lines)
print "  "+("OK" if resultWBWMQ else "FAIL")
good=good and resultWBWMQ
resultWBWMQr,failsWBWMQr=compare_runs(qpoly5files,qpoly4files,1,prec,lines)
print "  "+("OK" if resultWBWMQr else "FAIL")
goodr=goodr and resultWBWMQr

print" Comparing Polynomial versus Integration"
resultPI,failsPI=compare_runs(poly5files,int5files,prec,1,lines)
print "  "+("OK" if resultPI else "FAIL")
good=good and resultPI
resultPIr,failsPIr=compare_runs(poly5files,int5files,1,prec,lines)
print "  "+("OK" if resultPIr else "FAIL")
goodr=goodr and resultPIr

print "\nTest summary at absolute prec="+str(prec)+":"
print "         WB vs WM: "+("OK" if resultWBWM else "FAIL")
for s in failsWBWM:
    sys.stdout.write("           "+s)
print "   quad vs double: "+("OK" if resultQD else "FAIL")
for s in failsQD:
    sys.stdout.write("           "+s)
print "  WB vs WM (quad): "+("OK" if resultWBWMQ else "FAIL")
for s in failsWBWMQ:
    sys.stdout.write("           "+s)
print "      poly vs int: "+("OK" if resultPI else "FAIL")
for s in failsPI:
    sys.stdout.write("           "+s)
print "Overall absolute prec test: "+("PASSED" if good else "FAILED")


print "\nTest summary at relative prec="+str(prec)+":"
print "         WB vs WM: "+("OK" if resultWBWMr else "FAIL")
for s in failsWBWMr:
    sys.stdout.write("           "+s)
print "   quad vs double: "+("OK" if resultQDr else "FAIL")
for s in failsQDr:
    sys.stdout.write("           "+s)
print "  WB vs WM (quad): "+("OK" if resultWBWMQr else "FAIL")
for s in failsWBWMQr:
    sys.stdout.write("           "+s)
print "      poly vs int: "+("OK" if resultPIr else "FAIL")
for s in failsPIr:
    sys.stdout.write("           "+s)
print "Overall relative prec test: "+("PASSED" if goodr else "FAILED")


