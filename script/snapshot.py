#!/opt/local/bin/python
#use: module load other/SSSO_Ana-PyD/SApd_2.4.0_py2.7
import sys
import os
import numpy as np
import math
import subprocess
import argparse
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors
from matplotlib.backends.backend_pdf import PdfPages

nparmax=12

#functions:

#read step params from 
def get_step_pars(fname,stop_size=-1):
    fsize = os.stat(fname).st_size
    print( "fsize=",fsize)
    if(stop_size>=0):#sample from some specifed point in the middle of the file, instead of the end
        fsize=stop_size
    #print "fsize=",fsize
    bufsize = 15*(nparmax+5)
    with open(fname,'r') as f:
        #print "bufsize=",bufsize
        if bufsize > fsize:
            bufsize = fsize-1
        #print "bufsize=",bufsize
        #print "seeking=",fsize-bufsize
        data = []
        f.seek(fsize-bufsize)
        data.extend(f.readlines(bufsize))
        line=data[len(data)-3]
        print( "line is:",line )
        step= int(line.split()[0])
        pars= np.array(line.split()[5:])
        post= float(line.split()[1])
        #pars=pars[:-1]#last column is temp
    #print "step="+step
    return step,pars,post
# get_step_pars

#count the number of completed chains in the output file, based on number blank lines
def count_chains(fname):
    nblank=0
    chain_ends=[]
    pos=0
    with open(fname,'rb') as f:
        for line in f:
            if(line.startswith(b"\n")):
                #print "line:",line
                chain_ends.append(pos)
                #print "set: chain_ends["+str(int(nblank/2))+"]=",pos
                nblank+=1
            pos+=len(line)-1 #the "-1" is a hack because seek/tell/counting bytes all don't work together in an effective way. This does not seem to be related to the buffering issue on the tell() side with file iteration.
        chain_ends.append(pos)
        #bufsize = 15*(9+5)
        #print "testseeking=",chain_ends[0]-bufsize
        #data = []
        #f.seek(chain_ends[0]-bufsize)
        #data.extend(f.readlines(bufsize))        
        #print "testline is:",data[len(data)-3] 
    print( "chain_ends=",chain_ends)
    return nblank/2+1,chain_ends; #two blanks after completed chain.
#count_chains

#read arguments from outfile 
def get_flags(fname):
    d={}
    with open(fname) as f:
        line=f.readline()
        if(not line.count("flags=")):
           print( "I can't find the flags in the .out file: '",fname,"'")
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
                print( vals[0],"->",vals[1])
    return d
# get_flags
  
def get_param_text(basename):
    with open(basename+"_traj.dat") as file:
        parsline=file.readline()[1:]
        namesline=file.readline()[1:]
    namedpars=zip(namesline.split(),parsline.split())
    partext=""
    count=0;
    for pair in namedpars:
        count+=1
        if(count>5):
            count=1
            partext+="\n"
        partext+=pair[0]+"="+pair[1]+" "
        #print "count=",count,"pair=",pair
    return partext

def plot_lightcurve(basename,caption="",centerfrac=-1.0):
    #mimicing gnuplot: plot basename."_d_lcrv.dat" u 2:(-$3) ti "data" ,basename."_lcrv.dat" u 2:($4-$3) ti "" lw 2 w l,"" u 2:(-$4-$3) ti "model 1-sigma" lt 3 lw 2 w l, "" u 2:(-$3) ti "typical fit" lw 2 w l
    data=np.loadtxt(basename+"_d_lcrv.dat")
    model=np.loadtxt(basename+"_lcrv.dat")
    xd=data[:,1]
    xm=model[:,1]
    if(centerfrac>0):
        #find peak 90th percentile region
        percentile=.90
        y=np.append(data[:,2],model[:,2])
        y.sort()
        ycut=y[int(percentile*(len(y)-1))]
        xtop=np.array([l[1] for l in data if l[2]>ycut])
        xtop=np.append(xtop,np.array([l[1] for l in model if l[2]>ycut]))
        xcent=xtop.mean()
        xwid=xtop.std()
        print( "xcent=",xcent,"xwid=",xwid)
        xmin=xcent-centerfrac*xwid
        xmax= xcent+centerfrac*xwid
        print( "ycut=",ycut)
        print( xmin,"< x <",xmax)
        xd=np.ma.masked_where(np.any([xd < xmin , xd > xmax], axis=0), xd)
        xm=np.ma.masked_where(np.any([xm < xmin , xm > xmax], axis=0), xm)
    #plt.clf()
    fig=plt.figure()
    if(caption==""):ax1=fig.add_axes((.1,.1,.8,.8))
    else: ax1=fig.add_axes((.1,.2,.8,.7))
    ax1.plot(xd,-data[:,2],'+',color='darkorchid',label='data')
    ax1.plot(xm,-model[:,2]-model[:,3],'b-',label='_nolegend_')
    ax1.plot(xm,-model[:,2]+model[:,3],'b-',label='typical 1-sigma range')
    ax1.plot(xm,-model[:,2],'-',color='goldenrod',label='typical model')
    ax1.margins(y=0.05)
    ax1.legend()
    if(not caption==""):fig.text(.1,.05,caption)
                    
def plot_residual(basename,caption=""):
    #mimicing gnuplot: plot basename."_d_lcrv.dat" u 2:(-$3) ti "data" ,basename."_lcrv.dat" u 2:($4-$3) ti "" lw 2 w l,"" u 2:(-$4-$3) ti "model 1-sigma" lt 3 lw 2 w l, "" u 2:(-$3) ti "typical fit" lw 2 w l
    data=np.loadtxt(basename+"_d_lcrv.dat")
    xd=data[:,1]
    delta=-(data[:,2]-data[:,3])
    fig=plt.figure()
    if(caption==""):ax1=fig.add_axes((.1,.1,.8,.8))
    else: ax1=fig.add_axes((.1,.2,.8,.7))
    ax1.plot(xd,delta,'+',color='darkorchid',label='data')
    #ax1.errorbar(xd,delta,yerr=data[:,4],color='darkorchid',label='data')
    ax1.margins(y=0.05)
    ax1.legend()
    if(not caption==""):fig.text(.1,.05,caption)
                    
def plot_magmap(basename,caption="",var=""):
    #designed to mimic gnuplot:
    '''
    set view map
    set xtics border in scale 0,0 mirror norotate  offset character 0, 0, 0 autojustify
    set ytics border in scale 0,0 mirror norotate  offset character 0, 0, 0 autojustify
    set ztics border in scale 0,0 nomirror norotate  offset character 0, 0, 0 autojustify
    set rtics axis in scale 0,0 nomirror norotate  offset character 0, 0, 0 autojustify
    set cblabel "magnification" 
    set cbrange [*:10] noreverse nowriteback
    set palette rgbformulae -21, -22, -23
    set size square
    plot  basename."_mmap.dat" using 1:2:($3) with image, trajname u 3:4  lt 1 pointsize 0.1 lc rgb "black"
    rep ftrajname u 3:4  w l lc rgb "black" lw 0.1,trajname u 3:4 pointsize 0.02 lt 1 lc rgb "red"
    '''
    tol=1e-8;
    magmax=10.0;
    mdata=np.loadtxt(basename+var+"_mmap.dat")
    print( mdata.shape)
    x=np.unique(mdata[:,0])
    x=x[np.append([True],x[1:]-x[:-1]>tol)]
    dx=x[1:]-x[:-1]
    print( "dx limits=",max(dx),min(dx))
    y=np.unique(mdata[:,1])
    y=y[np.append([True],y[1:]-y[:-1]>tol)]
    dx=y[1:]-y[:-1]
    print( "dy limits=",max(dx),min(dx))
    nx=x.size
    ny=y.size
    print( "nx,ny,nx*ny,ndata",nx,ny,nx*ny,mdata[:,2].size)
    z=mdata[:,2]
    #z=(mdata[:,0]-min(x))/(max(x)-min(x))
    #z=(mdata[:,1]-min(y))/(max(y)-min(y))
    z=z.reshape(ny,nx)
    fig=plt.figure()
    if(caption==""):ax1=fig.add_axes((.1,.1,.8,.8))
    else: ax1=fig.add_axes((.1,.2,.8,.7))
    im=ax1.imshow(z, vmin=1.0,vmax=magmax,extent=(min(x),max(x),min(y),max(y)), origin='lower', cmap='hot_r')
    cb=plt.colorbar(im,ticks=np.arange(magmax)+1)
    #im=ax1.imshow(z, vmin=1.0,vmax=magmax,extent=(min(x),max(x),min(y),max(y)), origin='lower', cmap='hot_r')

    #superimpose the trajectory
    mtraj=np.loadtxt(basename+"_traj.dat")
    dtraj=np.loadtxt(basename+"_d_traj.dat")    
    #ax1.autoscale(False)
    ax1.autoscale(False)
    ax1.scatter( dtraj[:,2],dtraj[:,3],color='cyan',alpha=1,s=1.0,zorder=3)
    ax1.plot( mtraj[:,2],mtraj[:,3],'b-',scaley=False,scalex=False)
    if(not caption==""):fig.text(.1,.05,caption)
    
def make_plots(resultname, post="unknown"):
    pdf = PdfPages(resultname+'.pdf')
    caption="Log-posterior is "+str(post)+" for parameters:\n"+get_param_text(resultname)

    plot_lightcurve(resultname,caption)
    pdf.savefig()

    plot_residual(resultname,caption)
    pdf.savefig()

    plot_lightcurve(resultname,caption,1.5)
    pdf.savefig()

    plot_lightcurve(resultname,caption,0.5)
    pdf.savefig()

    plot_magmap(resultname,caption)
    pdf.savefig()

    plot_magmap(resultname,caption,"_z")
    pdf.savefig()

    plot_magmap(resultname,caption,"_zz")
    pdf.savefig()

    pdf.close()

##########
# setup stuff
# determine basename from first argument.

parser = argparse.ArgumentParser(description='Provide snapshot of chain state.')
parser.add_argument('fname', metavar='chain_file', type=str, 
                   help='chain file path')
parser.add_argument('datafile', metavar='data_file', type=str, nargs='?', default='',
                   help='lightcurve data file path')
parser.add_argument('-i','--ichain', metavar='chain_index', type=int, default=0,
                   help='index of chain to view')
parser.add_argument('-o','--offset', metavar='steps', action='store', type=int, default=0,
                   help='how far before end of chain to sample')
parser.add_argument('-e','--execname', metavar='executable_file', type=str, default="",
                   help='executable file path')
parser.add_argument('-v','--viewname', metavar='name', type=str, default="",
                   help='from a separate gleam-view run')
parser.add_argument('-p','--points', metavar='n', type=int, default=300,
                   help='specify number of sample points')

args = parser.parse_args()
print(args)

if(not args.viewname==""):
    make_plots(args.viewname)
    sys.exit()

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
print( "basename="+basename)

#get execname
#if(len(sys.argv)>2):
#    execname=sys.argv[2]
#else:
#    execname="../../src/gleam"
execname=""
if(len(args.execname)>0):execname=args.execname[0]
if(execname==""):execname=os.path.dirname(os.path.realpath(__file__))+"/../gleam"
print( "execname=",execname)

#get datafile
#if(len(sys.argv)>3):
#    datafile=sys.argv[3]
#else:
#    print "Need a method for automatically identifying the datafile, or provide it as third argument."
#    sys.exit()
datafile=args.datafile
print( "datafile=",datafile)

#get optional specification of which chain to take
#if(len(sys.argv)>4):
#    ichain=int(sys.argv[4])
#else:
#    ichain=-1

ichain=args.ichain
print( "ichain="+str(ichain))

npoints=args.points
print( "npoints=",npoints)

#get gnuplot script path
gpscript=sys.argv[0].replace(".py",".gp")

#determine arguments
d=get_flags(basename+".out")
eargs=""
for key in d.keys():
    eargs+=" -"+key+"="+d[key]
print( "eargs=",eargs)

#get pars and other info from .dat file
cc,cends=count_chains(fname)
print( "cends=",cends," cc=",cc)
iend=-1
if(ichain>=0 and cc>ichain):
    iend=cends[ichain]-args.offset
    cc=ichain
print( "iend=",iend," cc=",cc)

step,pars,post = get_step_pars(fname,iend)
resultname=basename+"_c"+str(cc)+"_"+str(int(math.floor(step/1000)))+"k"
parstr=""
for val in pars:
    parstr+=" "+str(val)
print(parstr)
parfile=resultname+".pars"
with open(parfile,"w") as f:f.write(parstr+"\n")

print( "step 0")

command= execname+" -view -mm_samples="+str(npoints)+" -stateFile="+parfile+" "+eargs+" "+datafile+" "+resultname
print( command)
print( command.split())
sys.stdout.flush()
subprocess.call(command.split())
print( "step 1")

make_plots(resultname,post)

command="gs "+resultname+".pdf"
print( command)
print( command.split())
sys.stdout.flush()
subprocess.call(command.split())
