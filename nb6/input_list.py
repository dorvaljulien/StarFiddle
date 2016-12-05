""" 
This file defines the "input_variables" function, which doesn't take 
any arguments and returns an I instance, which contains all default 
parameters for a Nbody6 run. These parameters can then be directly 
modified in I.
"""
import numpy as np
import os
#       Definition of input parameters, options & counters.
#       ---------------------------------------------------
## 
#        Input parameters
#        ****************
#

default_directory = os.path.join(os.environ["HOME"],"data/nbody6/test/")

def make_dic(**kwargs):
    return kwargs

#---------------------------------------------------------------------
def input_variables():
    class I:
        def __init__():
            pass



    I.executable = "nbody6"

#   Where the results will be placed
    I.directory = default_directory
    I.run_name="test"
    I.ForceErase=False
    I.path=I.directory+I.run_name
    I.Silent = False # print time steps while computing
    I.GPU=False # If the gpu version of nbody6, nbody6.gpu is available

    I.N       =  1000    # Number of objects (N_singles + 2*NBin)
    I.DELTAT  =  0.1     # Output time interval (N-body units).
    I.TCRIT   =  1.       # Termination time (N-body units).
    I.dt = I.DELTAT; I.tend = I.TCRIT
    

    I.model="King"
    # You can pick one of the following models:
    #         "King":         A King model
    #         "Plummer":      A Plummer sphere
    #         "Hubble":       Uniform sphere, with radial velocities: 
    #                         v=H.r
    #         "Uniform":      An uniform sphere, velocities are picked 
    #                         below escape velocities
    #         "Custom":       A snapshot of the initial state has to 
    #                         be provided: I.initial_state
    #         "NB6_uniform":  This is the built-in uniform sphere.
    #         "NB6_Plummer":  This is the built-in Plummer sphere.

    I.DontStart = False  # This will let the script prepare the directory
                         # with the raw input, input and fort.10 without
                         # actually starting the computation


    I.model_args = make_dic(N = 1000,
                            mass_range=[0.2,20], # Min and max masses.(Salpeter law)
                            virial_ratio = 0.5,  # On by default. Set to -1 to turn off.
                            # Hubble model
                            Hub = 1.16,           # For the model, vec(v)=Hub*vec(r)
                            virialise=False,      # True: The model will be virialised
                                                  # Hubble models are not virialised by default
                            # Plummer model
                            A0 = 0.1,             # The potential is  1/sqrt(A0**2 + r**2) 
                            truncature = 0.1,     # 0.4=> stars less likelky than 0.4% are rejected
                            #King model
                            W0=5
                            )



#   Custom model
    I.initial_state=""  # the initial state should be provided as a 
                        # snapshot without header: 
                        # m,x,y,z,vx,vz,vz    one line per star 

#   NB6_uniform.. Nothing capital to set

#   NB6_Plummer:
    I.AP      = -0.00     # Plummer scale factor 
                          # (N-body units; square saved in AP2).



#    NBODY6:

    I.KSTART =  1   # Control index 
                    # (1: new run; >1: restart; 3, 4, 5: new params).
    I.TCOMP   =  10000   # Maximum CPU time in minutes (saved in CPU).

#     INPUT

                          # Number of objects 
                          # (N_s + 2*N_b; singles + 3*NBIN0 < NMAX).
    I.NFIX    =  1        # Output frequency of data save or binaries 
                          # (options 3 & 6).
    I.NCRIT   =  25       # Final particle number (alternative 
                          # termination criterion).
    I.NRAND   =  10000    # Random number sequence skip.
    I.NBMAX  =  100       # Maximum number of neighbours (< LMAX - 5).
    I.NRUN    =  1        # Run identification index.

    I.ETAI    =  0.02     # Time-step parameter for irregular force 
                          # polynomial.
    I.ETAR    =  0.02     # Time-step parameter for regular force 
                          # polynomial.
    I.RS0     =  0.35     # Initial radius of neighbour sphere 
                          # (N-body units).
    I.DTADJ   =  0.20     # Time interval for parameter adjustment 
                          # (N-body units).
                          # Output time interval (N-body units).
                          # Termination time (N-body units).
    I.QE      =  8e-3     # Energy tolerance 
                          # (restart if DE/E > 5*QE & KZ(2) > 1).
    I.RBAR    =  5.5      # Virial cluster radius in pc 
                          # (set = 1 for isolated cluster).
    I.ZMBAR   =  0.5      # Mean mass in solar units 
                          # (=1.0 if 0; final depends on #20).


#--------------------------------------------------------------------------------------------
#------------------------        KZ VECTOR      ---------------------------------------------
#--------------------------------------------------------------------------------------------


    I.DTMIN   =  1.0E-04  # Time-step criterion for regularization search.
    I.RMIN    =  0.001    # Distance criterion for regularization search.
    I.ETAU    =  0.2      # Regularized time-step parameter (6.28/ETAU steps/orbit).
    I.ECLOSE  =  1.0      # Binding energy per unit mass for hard binary (positive).
    I.GMIN    =  1.0E-06  # Relative two-body perturbation for unperturbed motion.
    I.GMAX    =  0.001    # Secondary termination parameter for soft KS binaries.

# INPUT: if KZ(4) > 0
    I.DELTAS  = -0.00     # Output interval for binary search (in TCR; suppressed).
    I.ORBITS  = -0.00     # Minimum periods for binary output (level 1).
    I.GPRINT  = -0.00          # Perturbation thresholds for binary output (9 levels).

# DATA:

    I.ALPHAS  =  2.3      # Power-law index for initial mass function (used if #20 < 2).
    I.BODY1   =  10.0     # Maximum particle mass before scaling (KZ(20): solar mass).
    I.BODYN   =  0.2      # Minimum particle mass before scaling.
    NBIN0   =  0        # Number of primordial binaries (for IMF2 with KZ(20) > 1).
    I.NHI0    =  0        # Primordial hierarchies (may be needed in IMF if > 0).
    I.ZMET    =  0.002    # Metal abundance (in range 0.03 - 0.0001).
    I.EPOCH0  =  0.0      # Evolutionary epoch (in 10**6 yrs; NB! < 0 for PM evolution).
    I.DTPLOT  =  0.0      # Plotting interval for HRDIAG (N-body units; >= DELTAT).

# SETUP: if KZ(5) = 2

    I.APO     = -0.00     # Separation of two Plummer models (SEMI = APO/(1 + ECC).
    I.ECC     = -0.00          # Eccentricity of two-body orbit (ECC < 0.999).
    I.N2      = -0.00     # Membership of second Plummer model (N2 <= N).
    I.SCALE   = -0.00          # Second scale factor (>= 0.2 for limiting minimum size).

# If KZ(5) = 3

    I.APO3    = -0.00     # Separation between the perturber and Sun.
    I.ECC3    = -0.00          # Eccentricity of orbit (=1 for parabolic encounter).
    I.DMIN3   = -0.00     # Minimum distance of approach (pericentre).
    I.SCALE3  = -0.00          # Perturber mass scale factor (=1 for Msun).

# If KZ(5) = 4 

    I.SEMI4   = -0.00     # Semi-major axis (slightly modified; ignore if ECC > 1).
    I.ECC4    = -0.00          # Eccentricity (ECC > 1: NAME = 1 & 2 free-floating).
    I.M1_4    = -0.00     # Mass of first member (in units of mean mass).
    I.M2_4    = -0.00          # Mass of second member (rescaled total mass = 1).

# If KZ(5) >=6 and KZ(24) < 0

    I.ZMH     = -0.00     # Mass of single BH (in N-body units).
    I.RCUT    = -0.00          # Radial cutoff in Zhao cusp distribution (MNRAS, 278, 488).

# SCALE:

    I.NB6_virial =  0.5      # Virial ratio (Q = 0.5 for equilibrium).
    I.VXROT   =  0.0      # XY-velocity scaling factor (> 0 for solid-body rotation).
    I.VZROT   =  0.0      # Z-velocity scaling factor (not used if VXROT = 0).
    I.RTIDE   =  0.0      # Unscaled tidal radius (#14 >= 2; otherwise copied to RSPH2).
    I.SMAX    =  0.125    # Maximum time-step (factor of 2 commensurate with 1.0).

#     if (kz(24).gt.0)
#    M, X, V Initial subsystem (solar masses; membership = KZ(24)).
 
# XTRNL0:
# If KZ(14) = 2

    I.GMGPM   = -0.00     # Point-mass galaxy (solar masses, linearized circular orbit).
    I.RG0PM   = -0.00          # Central distance (in kpc).

# If KZ(14) = 3

    I.GMG     = -0.00     # Point-mass galaxy (solar masses).
    I.DISK    = -0.00          # Mass of Miyamoto disk (solar masses).
    I.A       = -0.00     # Softening length in Miyamoto potential (in kpc).
    I.B       = -0.00          # Vertical softening length (kpc).
    I.VCIRC   = -0.00     # Galactic circular velocity (km/sec) at RCIRC (=0: no halo).
    I.RCIRC   = -0.00          # Central distance for VCIRC with logarithmic potential (kpc).
    I.GMB     = -0.00     # Central bulge mass (Msun).
    I.AR      = -0.00          # Scale radius in gamma/eta model (kpc, Dehnen 1993).
    I.GAM     = -0.00     # gamma (in the range 0 =< gamma < 3).

    I.RG      = -0.00          # Initial X,Y,Z; DISK+VCIRC=0, VG(3)=0: A(1+E)=RG(1), E=RG(2).
    I.VG      = -0.00     # Initial cluster velocity vector (km/sec).

# If KZ(14)=3 or 4 

    I.MP      = -0.00          # Total mass of Plummer sphere (in scaled units).
                               # Plummer scale factor (N-body units; square saved in AP2).
    I.MPDOT   = -0.00          # Decay time for gas expulsion (MP = MP0/(1 + MPDOT*(T-TD)).
    I.TDELAY  = -0.00     # Delay time for starting gas expulsion (T > TDELAY).
    
# HOTSYS: If KZ(29) > 0

    I.SIGMA0  = -0.00          # Hot initial velocities in km/sec (CALL REFLCT suppressed).
    
# BINPOP: If KZ(8) = 1 or > 2

    I.NBIN0   =  0        # Number of primordial binaries (for IMF2 with KZ(20) > 1).
    I.SEMIB   =  0.005    # Max semi-major axis in model units (all equal if RANGE = 0).
    I.ECCB    =  -1.0     # Initial eccentricity (< 0 for thermal distribution).
    I.RATIOB  =  1.0      # Mass ratio M1/(M1 + M2); (= 1.0: M1 = M2 = <M>; not #20 > 1).
    I.RANGEB  =  5.0      # Range in SEMI for uniform logarithmic distribution (> 0).
    I.NSKIPB  =  5        # Binary frequency of mass spectrum (#20 < 2; body #1 first).
    I.IDORMB  =  0        # Indicator for dormant binaries (>0: merged components).
    I.ICIRCB  =  1        # Eigenevolution (=1: Kroupa & Mardling; =2: Kroupa 1995).
    I.MEANPB  =  0        # Added JDor - for lognormal period distribution - broken
    I.SIGMAPB =  0        # Added JDor - for lognormal period distribution - broken
    
# HIPOP: If KZ(8) > 0 and KZ(18) > 1

    I.SEMIH   = -0.00     # Max semi-major axis in model units (all equal if RANGE = 0).
    I.ECCH    = -0.00          # Initial eccentricity (< 0 for thermal distribution).
    I.RATIOH  = -0.00     # Mass ratio (= 1.0: M1 = M2; random in [0.5-0.9]).
    I.RANGEH  = -0.00          # Range in SEMI for uniform logarithmic distribution (> 0).
    I.ICIRCH  = -0.00     # Circularization & collision check (not implemented yet).
          
# INTIDE: If KZ(27) > 0  (currently suppressed)

    I.RSTAR   = -0.00          # Size of typical star in A.U.
    I.IMS     = -0.00     # idealized main-sequence stars.
    I.IEV     = -0.00          # idealized evolved stars.
    I.RMS     = -0.00     # Scale factor for main-sequence radii (>0: fudge factor).
    I.REV     = -0.00          # Scale factor for evolved radii (initial size RSTAR).
    
# INSTAR: if KZ(12) = 2  Input of stellar parameters on fort.12.
    
# CLOUD0: If KZ(13) > 0

    I.NCL     = -0.00     # Number of interstellar clouds.
    I.RB2     = -0.00          # Radius of cloud boundary in pc (square is saved).
    I.VCL     = -0.00     # Mean cloud velocity in km/sec.
    I.SIGMA   = -0.00          # Velocity dispersion (#13 > 1: Gaussian).
    I.CLM     = -0.00     # Individual cloud masses in solar masses (maximum MCL).
    I.RCL2    = -0.00          # Half-mass radii of clouds in pc (square is saved).
    

# CHAIN: If KZ(11) > 0 with ARchain

    I.CLIGHT  = -0.00     # Velocity of light in N-body units (e.g. 3.0D+05/VSTAR).
    I.NBH     = -0.00          # Number of BHs for special treatment (redundant but keep).
    I.IDIS    = -0.00     # Stellar disruption option (R_coll = (m/m_BH)^{1/3}*r^*).
    

    #---------------------------------------------------------------------


    #     Options KZ(J)
    #     *************
    I.KZ=np.zeros(51,dtype=np.int_)
#---------------------------------------------------------------------

    I.KZ[1] = 1  #  COMMON save unit 1 (=1: 'touch STOP'; =2: every 100*NMAX steps).

    I.KZ[2] = 1  #  COMMON save unit 2 (=1: at output; =2: restart if DE/E > 5*QE).

    I.KZ[3] = 0  #  Basic data unit 3 at output time (unformatted, frequency NFIX;
            #    =1/2: standard /and tail; 
            #    =3: tail only; 
            #    >3: cluster + tail).

    I.KZ[4] = 0  #  Binary diagnostics on unit 4 (# threshold levels = KZ(4) < 10).
            #  (currently suppressed in ksint.f.)

    I.KZ[5] = 8  #  Initial conditions (#22 =0; =0: uniform & isotropic sphere);
            #     =1: Plummer; 
            #     =2: two Plummer models in orbit, extra input;
            #     =3: massive perturber and planetesimal disk, extra input;
            #     =4: massive initial binary, extra input: A, E, M1, M2;
            #     =5: Jaffe model;
            #     =6: Zhao BH cusp model, extra input if #24 < 0: ZMH, RCUT.
            #     =7: King Model, see parameters at the beginning (HM JDor)
            #     =8: Uniform sphere, with Hubble-like velocity

    I.KZ[6] = 0  #  Significant & regularized binaries at main output (=1, 2, 3 & 4).

    I.KZ[7] = 1  #  Lagrangian radii 
            #     (>0: RSCALE; 
            #     =2, 3, 4: output units 6, 7, 12);
            #     >=2: half-mass radii of 50% mass, also 1% heavies, unit 6;
            #     >=2: Lagrangian radii for two mass groups on unit 31 & 32;
            #     >=2: geometric radii for three mass groups on unit 6;
            #     =5: density, rms velocity & mean mass on unit 26, 27 & 36;
            #     =6: pairwise values of mean mass and radii on unit 28.

    I.KZ[8] = 0  #  Primordial binaries 
            #     (=1 & >=3; >0: BINOUT; 
            #     >2: BINDAT; 
            #     >3: HIDAT;
            #     =4: Kroupa 1995 period distribution;
            #     >4: standard setup using RANGE & SEMI0).

    I.KZ[9] = 2  #  Individual bodies on unit 6 at main output (MIN(5**KZ9,NTOT)).

    I.KZ[10] = 0  #  Diagnostic KS output (>0: begin KS; >1: end; >=3: each step).

    I.KZ[11] = 0  #  Algorithmic Chain regularization and post-Newtonian (NBODY7).
             #  non-zero: PN for unpert KS or re-init ARchain (ksint.f);
             #     > 0: addition of initial BHs (binary/singles; scale.f);
             #     = -1: standard case of subsystem for ARchain (ksint.f);
             #     < -1: ARchain restricted to BH binary components (ksint.f).

    I.KZ[12] = 0  #  HR diagnostics of evolving stars (> 0; interval DTPLOT);
             #     =2: input of stellar parameters on fort.12 (routine INSTAR).

    I.KZ[13] = 0  #  Interstellar clouds (=1: constant velocity; >1: Gaussian).

    I.KZ[14] = 0  #  External force (
             #     =1: standard tidal field; 
             #     =2: point-mass galaxy;
             #     =3: point-mass + bulge + disk + halo + Plummer; 
             #     =4: Plummer).

    I.KZ[15] = 2  #  Triple, quad, chain (#30 > 0) or merger search (>1: more output).

    I.KZ[16] = 1  #  Updating of regularization parameters (>0: RMIN, DTMIN & ECLOSE);
             #     >1: RMIN expression based on core radius (experimental);
             #     >2: modify RMIN for GPERT > 0.05 or < 0.002 in chain.

    I.KZ[17] = 1  #  Modification of ETAI, ETAR (>=1) and ETAU (>1) by tolerance QE.

    I.KZ[18] = 1  #  Hierarchical systems (=1: diagnostics; =2: primordial; =3: both).

    I.KZ[19] = 0  #  Mass loss (
             #     =1: old supernova scheme; 
             #     =3: Eggleton, Tout & Hurley;
             #     >3: extra diagnostics).

    I.KZ[20] = 0  #  Initial mass function (
             #     =0: Salpeter type using ALPHAS; 
             #     =1: Scalo;
             #     =2, 4, 6: Kroupa; 
             #     =3, 5: Eggleton; 
             #     > 1: primordial binaries;
             #     =7: binary correlated m1/m2 also for brown dwarf IMF;
             #  Note: Use PARAMETER (MAXM=1) for setting BODY(1) = BODY10).

    I.KZ[21] = 1  #  Extra output (>0: MODEL #, TCOMP, DMIN, AMIN; >1: NESC by JACOBI).

    I.KZ[22] = 3  #  Initial m, r, v on #10 (
             #     =1: output; 
             #     >=2: input; 
             #     >2: no scaling;
             #     =2: m, r, v on #10 in any units; scaled to standard units;
             #                Note: choose #20 = 0 to avoid Salpeter IMF with scaling;
             #     =3: no scaling of input read on fort.10;
             #     =4: input from mcluster.c (no scaling; binaries if NBIN0 >0);
             #     =-1: astrophysical input (M_sun, km/s, pc) on unit #10).

    I.KZ[23] = 0  #  Escaper removal (
             #     >1: diagnostics in file ESC with V_inf in km/s);
             #     >=3: initialization & integration of tidal tail.

    I.KZ[24] = 0  #  Initial conditions for subsystem (M,X,V routine SCALE; KZ(24)= #);
             #  <0: ZMH & RCUT (N-body units) Zhao model (#5>=6).

    I.KZ[25] = 0  #  Velocity kicks for white dwarfs (
             #     =1: type 11 & 12; 
             #     >1: all WDs).

    # This option was suppressed
    # I.KZ[25] = 2  #  Partial reflection of KS binary orbit (GAMMA < GMIN; suppressed).

    I.KZ[26] = 0  #  Slow-down of two-body motion (>=1: KS; >=2: chain; =3: rectify).

    I.KZ[27] = 1  #  Tidal effects (
             #     =1: sequential; 
             #     =2: chaos; 
             #     =3: GR energy loss);
             #     =-1: collision detector, no coalescence, #13 < 0.

    I.KZ[28] = 0  #  GR radiation for NS & BH binaries (with #19 = 3; choice of #27);
             #  =4 and #27 = 3: neutron star capture (instar.f).

    I.KZ[29] = 0  #  Boundary reflection for hot system (suppressed).

    I.KZ[30] = 0  #  Multiple regularization (=1: all; >1: BEGIN/END; >2: each step);
             #  =-1: CHAIN only; =-2: TRIPLE & QUAD only. 

    I.KZ[31] = 1  #  Centre of mass correction after energy check.

    I.KZ[32] = 0  #  Increase output intervals & SMAX based on single particle energy.

    I.KZ[33] = 0  #  Histograms at main output (>=1: STEP; =2: STEPR, NBHIST & BINARY).

    I.KZ[34] = 0  #  Roche-lobe overflow (=1: Roche & Synch; =2: Roche & BSE synch).

    I.KZ[35] = 0  #  Time offset (global time from TTOT = TIME + TOFF; offset = 100).

    I.KZ[36] = 0  #  Step reduction for hierarchical systems (suppressed).

    I.KZ[37] = 1  #  Neighbour additions in CHECKL (>0: high-velocity; >1: all types).

    I.KZ[38] = 0  #  Force polynomial corrections (=0: standard, no corrections;
                  #     =1: all gains & losses included;
                  #     =2: small FREG change skipped;
                  #     =3: fast neighbour loss only).

    I.KZ[39] = 0  #  No unique density centre (skips velocity modification of RS(I)).

    I.KZ[40] = 2  #  Neighbour number control (=1: increase if <NNB>  <  NNBMAX/2);
             #  >=2: fine-tuning at NNBMAX/5; =3: reduction of NNBMAX.

    I.KZ[41] = 0  #  Pre-mainsequence stellar evolution (only solar metallicity).

    I.KZ[42] = 0  #  Kozai diagnostics on fort.42 (=1: frequency 100 & EMAX > 0.99).

    I.KZ[43] = 0  #  Small velocity kick after GR coalescence (=1, =3; NBODY7 only),
             #  =2: BH accretion of disrupted star, KSTAR >= 10.

    I.KZ[44] = 0  #  Plotting file for main cluster parameters on fort.56 (OUTPUT).

    I.KZ[45] = 0  #  Plotting file for BH (NAME = 1 or 2) on unit 45 (routine BHPLOT);
             #  primordial BH defined by INSTAR; membership = KZ(24).

    I.KZ[46] = 0  #  Reserved for data analysis project on NBODY6++.

    I.KZ[47] = 0  #  Reserved for data analysis project on NBODY6++.

    I.KZ[48] = 0  #  GPU initialization of neighbour lists and forces (FPOLY0).

    I.KZ[49] = 0  #  Post-Newtonian perturbations included in KS (dir Block).
    # ---------------------------------------------------------------------


    

    # Path of the working directory with the nbody6 executable
    I.nb6_directory="/home/dorval/nbody6/main_install/work/"


    return I

 
 #  ---------------------------------------------------------------------

 #  Output counters
 #  ***************

 #  ---------------------------------------------------------------------
 #  NSTEPI  Irregular integration steps.
 #  NSTEPR  Regular integration steps.
 #  NSTEPU  Regularized integration steps.
 #  NNPRED  Coordinate & velocity predictions of all particles.
 #  NBPRED  Coordinate & velocity prediction of neighbours (NNB counted).
 #  NBCORR  Force polynomial corrections.
 #  NBFULL  Too many neighbours with standard criterion.
 #  NBVOID  No neighbours inside 1.26 times the basic sphere radius.
 #  NICONV  Irregular step reduction (force convergence test).
 #  NBSMIN  Retained neighbours inside 2*RS (STEP < SMIN).
 #  NLSMIN  Small step neighbours selected from other neighbour lists.
 #  NBDIS   Second component of recent KS pair added as neighbour (#18).
 #  NBDIS2  Second component of old KS pair added as neighbour (#18 > 1).
 #  NCMDER  C.m. values for force derivatives of KS component.
 #  NBDER   Large F3DOT corrections not included in D3 & D3R.
 #  NFAST   Fast particles included in LISTV (option 37).
 #  NBFAST  Fast particles included in neighbour list (option 37).
 #  NBLOCK  Number of blocks (block-step version).
 #  NBLCKR  Number of regular force blocks.
 #  NMDOT   Mass loss events (option 19).
 #  NBSTAT  Diagnostic data on binary interactions (option 4).
 #  NKSTRY  Two-body regularization attempts.
 #  NKSREG  Total KS regularizations.
 #  NEWKS   Enforced KS regularization using wider criterion (~8 > 0).
 #  NKSHYP  Hyperbolic KS regularizations.
 #  NKSPER  Unperturbed KS binary orbits.
 #  NPRECT  Initialization of NKSPER after exceeding 2*10**9.
 #  NKSREF  Partial reflections of KS binary (option 25; suppressed).
 #  NKSMOD  Slow KS motion restarts (option 26).
 #  NTTRY   Search for triple, quad & chain regularization or mergers.
 #  NTRIP   Three-body regularizations (option 15).
 #  NQUAD   Four-body regularizations (option 15).
 #  NCHAIN  Chain regularizations (options 15 & 30).
 #  NMERG   Mergers of stable triples or quadruples (option 15).
 #  NEWHI   New hierarchical systems (counted by routine HIARCH).
 #  NSTEPT  Triple regularization integration steps (option 15).
 #  NSTEPQ  Quadruple regularization integration steps (option 15).
 #  NSTEPC  Chain regularization steps (# DIFSY calls).
 #  NDISS   Tidal dissipations at pericentre (option 27).
 #  NTIDE   Tidal captures from hyperbolic motion (option 27).
 #  NSYNC   Number of synchronous binaries (option 27).
 #  NCOLL   Stellar collisions (option 27).
 #  NSESC   Escaped single particles (option 23).
 #  NBESC   Escaped binaries (option 23).
 #  NMESC   Escaped mergers (options 15 & 23).
 #  NRG     Red giants.
 #  NHE     Helium stars.
 #  NRS     Red supergiants.
 #  NNH     Naked Helium stars.
 #  NWD     White dwarfs.
 #  NSN     Neutron stars.
 #  NBH     Black holes.
 #  NBS     Blue stragglers.
 #  ---------------------------------------------------------------------


 #  Stellar evolution types
 #  ***********************

 #  ---------------------------------------------------------------------
 #  0       Low main sequence (M < 0.7).
 #  1       Main sequence.
 #  2       Hertzsprung gap (HG).
 #  3       Red giant.
 #  4       Core Helium burning.
 #  5       First AGB.
 #  6       Second AGB.
 #  7       Helium main sequence.
 #  8       Helium HG.
 #  9       Helium GB.
 # 10       Helium white dwarf.
 # 11       Carbon-Oxygen white dwarf.
 # 12       Oxygen-Neon white dwarf.
 # 13       Neutron star.
 # 14       Black hole.
 # 15       Massless supernova remnant.
 # 19       Circularizing binary (c.m. value).
 # 20       Circularized binary.
 # 21       First Roche stage (inactive).
 # 22       Second Roche stage.
 #  ---------------------------------------------------------------------


