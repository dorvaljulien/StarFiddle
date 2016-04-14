"""
This module allows a user to generate an input file for NBODY6 
and launch the integration in an user-friendly manner.
"""

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import cPickle as pk
import collections
import os
import subprocess
import datetime
import inspect

from readcol import *
import amuse.lab as AL
import amuse.units.units as U
import amuse.units.nbody_system as nbU

from nb6.input_list import input_variables
from nb6.pysnap import PySnap
import animation3d as anim3d
import anim_binary_histogram as animbin
import anim_binary_parameters as animparam
import cluster_models as CM
from nb6.miscnb6 import read_run, snapname
from distutils.util import strtobool

import binaries

AU_to_parsec=4.84e-6 

#---------------------------------------------------------------------
#--------------------------------------------------------------------
#--------------------------------------------------------------------
#--------------------------------------------------------------------
#--------------------------------------------------------------------
def user_yes_no_query(question):
    sys.stdout.write('%s [y/n]\n' % question)
    while True:
        try:
            return bool(strtobool(raw_input().lower()))
        except ValueError:
            sys.stdout.write('Please respond with \'y\' or \'n\'.\n')


def LaunchNbody6(**kwargs ):
    """
    Initiates a Nbody6 run. All arguments in input_list.py can be set.
    Return a container of all arguments.
    """
    I=input_variables()
    for variable, value in kwargs.iteritems():
        if hasattr(I,variable):
            setattr(I,variable,value)
        if I.model_args.has_key(variable):
            I.model_args[variable]=value
    I = launch_from(I)
    return I
        



def launch_from(I):
    """
    Take a container of arguments for a NBODY6 run, prepare the directory, launch
    NBODY6 then watch the evolution.
    """
    # Adding final "/"
    if I.directory[-1]!="/":
        I.directory=I.directory+"/"
    I.path=I.directory+I.run_name
    if I.path[-1]!="/":
        I.path=I.path+"/"

    # Dealing with some stuff
    try:
        I.TCRIT=I.tend
    except AttributeError:
        pass
    try:
        I.DELTAT=I.dt
    except AttributeError:
        pass
    

    # We first create the target directory, or clean it if it already exists
    if os.path.isdir(I.path):
        if not I.ForceErase:
            reply = user_yes_no_query("A run already exists at:\n   "+I.path+
                                      "\nDo you want to erase it ?")
        if I.ForceErase or reply:
            subprocess.call("clean_output.sh "+I.path,shell=True) 
        else:
            print "Exiting..."
            sys.exit(1)
    else:
        subprocess.call("mkdir "+I.path+" 2>/dev/null",shell=True)
    

    def generate_King():
        if not I.Silent: 
            print "Generating King model..."
        return CM.King( **I.model_args )
    def generate_Plummer():
        if not I.Silent: 
            print "Generating Plummer model..."
        return CM.Plummer( **I.model_args )
    def generate_Hubble():
        if not I.Silent: 
            print "Generating Hubble velocities model..."
        return CM.Hubble( **I.model_args )
    def generate_Uniform():
        if not I.Silent: 
            print "Generating Uniform sphere model... "
        return CM.Uniform( **I.model_args )
    def generate_Custom():
        if not I.Silent: 
            print "Loading custom model..."
        subprocess.call("cp "+I.initial_state+" "+
                        I.path+"fort.10",shell=True)
        with open(I.initial_state) as f:
            I.N=len(f.readlines())


    # We get the appropriate initial model:
    model_generator={"King":generate_King,
                    "Plummer":generate_Plummer,
                    "Hubble":generate_Hubble,
                    "Uniform":generate_Uniform,
                    "Custom":generate_Custom
                    # "NB6_Uniform":generate_NB6_Uniform,
                    # "NB6_Plummer":generate_NB6_Plummer
                 }

    if I.model is "Custom" or isinstance(I.initial_state, PySnap ):
        if not I.Silent:
            print "Loading custom model..."
        if isinstance(I.initial_state, PySnap ):
            I.initial_state = I.initial_state.WriteSnapshot(NB6_InputFormat=True)
        subprocess.call("cp "+I.initial_state+" "+
                        I.path+"fort.10",shell=True)
        with open(I.initial_state) as f:
            I.N=len(f.readlines())
    else:
        try:
            S = model_generator[I.model]() # We call the appropriate function
            if I.DontSplitBinaries:
                print "Dumping binaries population in "+I.path+ " as BinPop.pk"
                pk.dump(S.BinPop, open(I.path+"BinPop.pk","w") )
            if hasattr(S,"scaling"):
                with open(I.path+"scaling", "r") as f:
                    f.write(S.scaling)
            S.WriteSnapshot(I.path+"fort.10",NB6_InputFormat=True)

        except KeyError:
            raise Exception("Please pick an initial model:\n"
                            "Hubble, King, Plummer, Uniform, Custom")

    # Tell nbody6 to read the model on fort.10
    I.KZ[22]=3
        
    # We now write the raw_input file which is readable by nbody6:
    raw_filename=I.path+"raw_input"
    f=open(raw_filename,"w")

    # ----------------------------------------------------------------
    # ----------------------------------------------------------------
    # ----------------------------------------------------------------
    # Now for the fun part, let's write the raw_input file
    

    #NBODY6
    list1=[str(I.KSTART),str(I.TCOMP)]
    f.write(" ".join(list1)+"\n")
    
    #INPUT
    list2=[str(I.N),str(I.NFIX),str(I.NCRIT),
           str(I.NRAND),str(I.NBMAX),str(I.NRUN)]
    f.write(" ".join(list2)+"\n")
    
    list3=[str(I.ETAI),str(I.ETAR),str(I.RS0),
           str(I.DTADJ),str(I.DELTAT),str(I.TCRIT),
           str(I.QE),str(I.RBAR),str(I.ZMBAR)]
    f.write(" ".join(list3)+"\n")
    
    listk1=[str(I.KZ[i+1]) for i in range(10)]
    listk2=[str(I.KZ[i+1]) for i in range(10,20)]
    listk3=[str(I.KZ[i+1]) for i in range(20,30)]
    listk4=[str(I.KZ[i+1]) for i in range(30,40)]
    listk5=[str(I.KZ[i+1]) for i in range(40,50)]
    

    f.write(" ".join(listk1)+"\n")
    f.write(" ".join(listk2)+"\n")
    f.write(" ".join(listk3)+"\n")
    f.write(" ".join(listk4)+"\n")
    f.write(" ".join(listk5)+"\n")
    
    list4=[str(I.DTMIN),str(I.RMIN),str(I.ETAU),
           str(I.ECLOSE),str(I.GMIN),str(I.GMAX)]
    f.write(" ".join(list4)+"\n")
    
    
    if I.KZ[4] > 0:
        list5=[str(I.DELTAS),str(I.ORBITS),str(I.GPRINT)]
        f.write(" ".join(list5)+"\n")

    #DATA
    list6=[str(I.ALPHAS),str(I.BODY1),str(I.BODYN),str(I.NBIN0),
           str(I.NHI0),str(I.ZMET),str(I.EPOCH0),str(I.DTPLOT)]
    f.write(" ".join(list6)+"\n")


    #SETUP
    if  I.KZ[5]==2:
        list7=[str(I.APO),str(I.ECC),str(I.N2),str(I.SCALE)]
        f.write(" ".join(list7)+"\n")

    if I.KZ[5]==3:
        list8=[str(I.APO3),str(I.ECC3),str(I.DMIN3),str(I.SCALE3)]
        f.write(" ".join(list8)+"\n")

    if I.KZ[5]==4:
        list9=[str(I.SEMI4),str(I.ECC4),str(I.M1_4),str(I.M2_4)]
        f.write(" ".join(list9)+"\n")

    if I.KZ[5]>=6 and I.KZ[24]<0:
        list10=[str(I.ZMH),str(I.RCUT)]
        f.write(" ".join(list10)+"\n")


    # SCALE
    list11=[str(I.NB6_virial),str(I.VXROT),str(I.VZROT),
            str(I.RTIDE),str(I.SMAX)]
    f.write(" ".join(list11)+"\n")


    # XTRNL0
    if I.KZ[14]==2:
        list12=[str(I.GMGPM),str(I.RG0PM)]
        f.write( " ".join(list12)+"\n")

    if I.KZ[14]==3:
        list13=[str(I.GMG),str(I.DISK),str(I.A),str(I.B),
                str(I.VCIRC),str(I.RCIRC),str(I.GMB),
                str(I.AR),str(I.GAM)]
        f.write(" ".join(list13)+"\n")

        list14=[str(I.RG),str(I.VG)]
        f.write(" ".join(list14)+"\n")

    if I.KZ[14]==3 or I.KZ[14]==4:
        list15=[str(I.MP),str(I.AP),str(I.MPDOT),str(I.TDELAY)]
        f.write(" ".join(list15)+"\n")

    # HOTSYS
    if I.KZ[29]>0:
        list16=[str(I.SIGMA0)]
        f.write(list16+"\n")

    # BINPOP
    if I.KZ[8]==1 or I.KZ[8]>2:
        list17=[str(I.SEMIB),str(I.ECCB),str(I.RATIOB),
                str(I.RANGEB),str(I.NSKIPB),str(I.IDORMB),
                str(I.ICIRCB),str(I.MEANPB),str(I.SIGMAPB)]
        f.write(" ".join(list17)+"\n")

    # HIPOP
    if I.KZ[8]>0 and I.KZ[18]>1:
        list18=[str(I.SEMIH),str(I.ECCH),str(I.RATIOH),
                str(I.RANGEH),str(I.ICIRCH)]
        f.write(" ".join(list18)+"\n")

    # INTIDE
    if I.KZ[27]>0:
        list19=[str(I.RSTAR),str(I.IMS),str(I.IEV),
                str(I.RMS),str(I.REV)]
        f.write(" ".join(list19)+"\n")

    # CLOUD0
    if I.KZ[13]>0:
        list20=[str(I.NCL),str(I.RB2),str(I.VCL),str(I.SIGMA),
                str(I.CLM),str(I.RCL2)]
        f.write(" ".join(list20)+"\n")

    # CHAIN
    if I.KZ[11]>0:
        list21=[str(I.CLIGHT),str(I.NBH),str(I.IDIS)]
        f.write(" ".join(list21)+"\n")

    datestring = datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")
    
    f.write("\n")
    f.write("\n")

    f.write("Automatically generated for nbody6, "+datestring+
            "\nJulien Dorval - dorvaljulien@gmail.com")

    f.close()

    # ----------------------------------------------------------------
    # ----------------------------------------------------------------
    # ----------------------------------------------------------------

    if not I.DontStart:
        # We now call for nbody6, then keep watching the log file
        cwd=os.getcwd()
        os.chdir(I.path)
        if I.GPU:
            exec_name="nbody6.gpu"
        else:
            exec_name=I.executable
            
        if I.Silent:
            subprocess.call(exec_name+" <"+raw_filename+
                            " 2>1 1>log", shell = True)
        else:
            subprocess.call(exec_name+" <"+raw_filename+
                            " 2>&1 | tee "+I.path+"log"+
                            " | grep \" T =\"",
                            shell=True)   

        # We dump the input as a text file:
        with open(I.path+"input","w") as f:
            for variable in vars(I):
                f.write(variable+"\t\t"+str(getattr(I,variable))+"\n")

        # We go back to the working directory
        os.chdir(cwd)

    return I



