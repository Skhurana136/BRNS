# _________________________________________________________________________

# AUTOMATIC CODE GENERATOR (ACG) FOR CONSTRUCTING USER-DEFINED
# BIOGEOCHEMICAL REACTION NETWORKS
# 
# version 1.2
# COPYRIGHT (c) 2001 P.A.G. Regnier 
# All Rights Reserved
# 
# Research Unit on Biogeochemical Systems Dynamics
# Department of Geochemistry, Utrecht University, 
# The Netherlands
# __________________________________________________________________________
# __________________________________________________________________________
# INPUT TYPES
# 
# OOO : Sections that should be modified by the user
# OOO : Sections that should NOT be modified by the user

# OOO : comments 
# OOO : Maple input
# OOO : Maple output (appears only after you have executed the spreadsheet) 
# OOO  : Maple input entries that have to be specified by the user 
# 
# WWW  : Hyperlink to the Knowledge Book
# 
# _________________________________________________________________________ 

#  Maple specific info
restart ;
#with(Spread) :
 precision := double :
#  Summary and governing equations
#  Governing equation 
# The governing equation solved is of the form :
# theta(x)*dC/dt = -d(Fdiff+Fadv)/dx + R,
# where 
# - dX/dp is the partial derivative of X with respect to p [M/L^3/[p]]
# - theta(x) = A(x)*por(x) for dissolved, and A(x)*(1-por(x)) for solid species [L^2]
# - Fdiff = -(D(x)*theta(x)*dC/dx) [M/L^2/T]
# - Fadv  = (v*theta(x)*C)[M/L^2/T]
# - A is the cross section area [L^2], por is porosity [-]
# - v is the advection velocity [L/T], the sum of 
# the water flow velocity (vwat = water flow q in [L^3/T]/theta(x), solutes only) and 
# an advective velocity w acting upon both solids and solutes (e.g. movement due to fixed reference frame)
# - D is the effective diffusion coefficient [L^2/T], the sum of 
# the molecular diffusion coefficient (Dmol, solutes only), 
# the bioturbation coefficient (Db) and 
# the dispersion coefficient (Disp = aL*|vwat|, solutes only, aL in [L]). 
# - R is the sum of the reaction terms [M/L^3/T] as specified in the reaction network established below
#  Caveats 
# - for real numbers you should add a point after the number.
# - make sure all units match. 
# e.g. rates and dXdt, or flux boundary conditions and concentrations
# for the conversion of solid to solute units, 
# one may define temporary variables that can be used in the rate laws below:
# s_dens := 2.5; # solid density in [g/cm_solid^3]
# sd := 1000. * s_dens * (1. - por(j)) / por(j); 
# the factor 1.d03 converts cm^3 to liter, e.g. [g/cm^3] to [g/l]. 
# note that you need to refer to porosity exactly as por(j).
#  Physics - Parameters
# The list of length nparphys & nparphys2 is given by phys_name & phys_name2; the values collected in phys_val & phys_val2.
#  Spatial and temporal domain size
# tot_time: total length of simulation (T, e.g. years)
# depth_max: spatial extent, e.g. total depth of simulation (L, e.g. cm)
tot_time := 5000. ;
depth_max := 100.0 ;
#  Transport coefficients
# vaL: longitudinal dispersivity (L). The dispersion coefficient is calculated as Disp = aL*|vwat| 
# viq: depth dependency of flow. 0: constant, else changing with depth (i.e. distance)
# vq0: water flow (either constant or else at x=0, L^3/T)
# viw: depth dependency of w. 0: constant, else changing with depth (i.e. distance)
# vw0: advection velocity working on solids and solutes (either constant or else at x=0, L/T)
# viDb: depth dependency of Db. 0: constant, else changing with depth (i.e. distance)
# vDb0: bioturbation coefficient working on solids and solutes (either constant or else at x=0, L^2/T)
# 
# - user specified profiles of q(x), w(x) or Db(x) must be defined in the fortran subroutines advcoeff.f and diffcoeff.f, respectively.
# - the molecular diffusion coefficients for the species used are specified further below
val := 5.0 ;
viq := 0 ;
vq0 := 0. ;
viw := 0 ;
vw0 := 100.0 ;
viDb := 0 ;
vDb0 := 0.0 ;
#  Porosity profile and cross section area
# vipor: depth dependency of porosity. 0: constant, else changing with depth (i.e. distance)
# vpor0: porosity value (either constant or else at x=0)
# viarea: depth dependency of cross section area. 0: constant, else changing with depth (i.e. distance)
# varea0: cross section area (either constant or else at x=0, L^2)
# - user specified profiles of por(x) or area(x) must be defined in the fortran subroutine porarea.f.
vipor := 0 ;
vpor0 := 0.4 ;
viarea := 0 ;
varea0 := 1.0 ;
#  Temperature and salinity
# T_C: temperature (Celsius)
# S : salinity (PSU)
# T and S are used to calculate the molecular diffusion coefficients
# T: absolute temperature (Kelvin)
T_C := 20. ;
S := 35.00 ;
T := T_C + 273.15 ;
#  Grid and discretization
# Dt: time step of numerical integration (T). 
# note that there are options defined in drivervalues.f that allow automatical selection of timestep
# nnodes: number of nodes in the spatial domain. 
# for a regular grid, the grid spacing of the concentration profile is then depth_max/(nnodes-1)
# vigridv: type of grid. 
# 0 is for regular, evenly spaced grid, 
# else the user needs to specify the grid in the fortran subroutine gridsetup.f
Dt := 0.005 ;
nnodes := 51;
vigrid := 0 ;
#  Reaction Network - Size and Variables
#  Size of reaction network
# nsolids : number of solid species 
# ndissolved : number of dissolved species 
# ncompo : total number of species 
# nreactions : total number of reactions (including equilibrium rxns) 
# neqrxns : number of equilibrium reaction
nsolids := 2 ;
ndissolved := 4 ;
ncompo := nsolids + ndissolved ;
nreactions := 5 ;
neqrxns := 0 ;
#  List of variables
# variables: list of variables to model. example: 
# listsolids: species number which is a SOLID species. 
# note: all other variables are temporary and are NOT parsed to the ACG
# 
# Example: 
# variables:=[O2, so4, MnOx, FeOx, hco3, co3, hplus, hs];
# listsolids: = [3,4];
 variables := [ch2o,c_bioav,baco,bacs,o2,so4] ;
listsolids := [3,4] ;
#  Biogeochemistry - Rate laws
# Definition of kinetic rate laws
# rate.i : array of rates. 
# - For equilibrium rate expression, a kinetic rate MUST be specified as well. It will be overwritten in the equilibrium section below, but you need it as space holder and stoichiometry. Furthermore, the steadystate module uses detailed balancing method with fast kinetics. Therefore, in the example below, kf will have to be defined as a rather large number and kb = kf*Keq
# note: all other variables are temporary and are NOT parsed to the ACG
# - conditional statements: if a rate law depends on a conditional statement you need to make use of the subroutine switches.f. Example: dissolution (Rd) is only to take place at undersaturation, thus  Rd= f(saturation). If saturation > 1, Rd>0, else Rd=0. This canbe implemented as Rd:=k*H1*saturation, where H1 is toggled between 0 and 1. Rather than giving the condition here in maple, for now you need to do this in "switches.f", where you program the conditions for H1, e.g. H1=0, If (A*B/K>1) then H1 = 1.
#  
# example: 
# rate1 := 1000.*O2*hs; # rate law for 2O2 + HS -> SO4 + Hplus
# rate2 := kf*hplus*co3 - kb*hco3; # kinetic rate law for HCO3 = CO3 + Hplus (equilibrium)
# 
#  Primary redox reactions WWW

#  Sulfide Re-oxidation reactions
rate1:= muc*(ch2o-c_bioav) ;
rate2:= muo*baco*c_bioav/(kmoc+c_bioav)*o2/(kmo+o2) ; #*(1-(baco+bacn+bacs)/baccmax);
rate3:= deco*baco*(1-bacmin/baco) ;
rate4:= mus*bacs*c_bioav/(kmsc+c_bioav)*so4/(kms+so4)*kmoinh/(kmoinh+o2); #*(1-(baco+bacn+bacs)/baccmax); 
rate5:= decs*bacs*(1-bacmin/bacs) ;

#  Biogeochemistry - Stoichiometry

# Stoichiometry of the biogeochemical reactions
# d.sp.dt : rates of change of sp due to the sum of biogeochemical reactions
# note that rateX must be referred to as rX  
#  
# example:
# dO2dt := -2*r1;
# dhco3dt = -r2; 
s_dens := 2.5 ;
#SD := 1.0e3 * s_dens * (1.0 - por(j)) / por(j);
# 
dch2odt := -r1 ;
dc_bioavdt := r1 - r2/yieldo*SD - r4/yields*SD;
dbacodt := r2 - r3 ;
dbacsdt := r4 - r5 ;
do2dt := -r2/yieldo*(1-yieldo)*SD ;
dso4dt := -r4/yields*(1-yields)*0.5*SD ;


#  Biogeochemistry - Equilibria
# Specification of equilibrium constraints
# eqrxnsId : set of kinetic reactions which are overuled by a thermodynamic constraint 
# equilibriumseqns[i] : Equilibrium constraint for reaction i  
# 
# example:
# eqrxnID := [r2,rX];
# equilibriumeqns[1] := hplus*co3 - Keq*hco3;
# equilibriumeqns[2] := ...;
eqrxnId := [] ;
#equilibriumeqns[1] := b - k3*c ;
#  Biogeochemistry - Parameters
# Values of rates constants and parameters 
# In this section, all parameters defined in section 'Rate laws' should be defined.
# nparam: number of parameters to define 
# The list is given by bio_name; the values collected in bio_val.
# note that for double precision, 10 should be written as 10.
# 
# example:
# nparam:=4;
# bio_name:=[kmo2hs,kf,kb,Keq];
# vkf :=1.0*10^(5);
# vKeq:=1.0*10^(-10.4);
# vkb :=vkf*vKeq;
# bio_val:=[1000.,vkf,vkb,vKeq];
# 
# 
nparam := 14 ;
bio_name := [muc,bacmin,   kmoc,muo,deco,kmo,yieldo,   mus,decs,yields,kms,   kmsc,kmoinh,SD] ; 
bio_val :=[2.0,1.0e-9,   1.0e-6,6.0,0.01,1.0e-6,0.2,   0.6,0.01,0.02,1.0e-6,    kmoc,kmo,s_dens*(1-vpor0)/vpor0*1000] ; 

# Switches
# Switches can be used in the rate equations. Specify in nswitches, how many switches are in use, name them and define the switch expressions. The switch names must also appear in bio_name and must be assigned a dummy value there. The switch equals 1 if the switch expression is >0, 0 otherwise. To reference the coordinates in the domain, use x_pos, y_pos and z_pos.
nswitches := 0 ;
switchlist := [sw1,sw2,sw3] ;
switchcrit := [(bio*dissb-0.25),dissc,dissc] ;
# 
#  Transport - Molecular Diffusion
# Spec ification of the molecular diffusion coefficients
# diffdata: molecular diffusion coefficient at 0 degree celsius (cm^2/yr)
# alphadata: temperature dependence of the diffusion coefficient (1/K) 
# the in situ molecular diffusion coefficient, corrected for tortuosity is calculated as:
# D(T,sal) = [(0.95-0.001*sal)* D(T=0,sal=0)*(1 + muc*T[C])]/(1-ln(por^2))  
# example:
# diffdata := [100., 304.,0.,0.,100.];
# alphadata:= [0.006, 0.04,0.,0.,.0.05];
diffdata := [1.0,1.0,0.0,0.0,   1.0,1.0] ;
alphadata := [0.06,0.06,0.06,0.06,0.06,0.06] ; 

#  Transport - Boundary Conditions
# Specification of upper and lower boundary conditions for each species.
# There are 3 options 
# 0. kmnown concentration (Dirichlet, M/L^3)
# 1. kmnown concentration gradient (Neumann, M/L^3*L)
# 2. kmnown total (diffusive and advective) flux (Robin, M/L^2/T)
# technical note: option 1 and 2 involve ghost points outside the domain. If the mixing parameters vary with depth, one needs to assign a mixing intensities at the ghost points. By default this is done by linear extrapolation. To ovrwrite this, the user has to edit gridsetup.f and advdiffcoeff.f (both at the bottom; explanations are given there)
# type_up: array defining the type of condition for each species at the upper boundary (0, 1 or 2) 
# bnddata_up: array containing the values specified at the upper boundary.
# type_down: array defining the type of condition for each species at the lower boundary 
# bnddata_down: array containing the values specified at the upper boundary.
# 
# example:
# type_up := [0,0,2];
# bnddata_up := [1.4,0.001,0.001];
# type_down := [1,1,1];
# bnddata_down := [0.,0.,0.];
type_up := [2,2,2,2,2,2] ;
bnddata_up := [0.75e-3*40.0,0.0,0.0,0.0,0.25e-3*40.0,0.25e-3*40.0] ; 
type_down := [1,1,1,1,1,1] ;
bnddata_down := [0.0,0.0,0.0,0.0,0.0,0.0] ;
#  Transport - Identify species that are not transported
# nrnotransp: number of species that are not getting transported at all
# notransp: list of species that are not getting transported at all 
# (i.e. for which the above defined dispersion and advection velocities are not used!)
# 
# example:
# nrnotransp:=3;
# notransp := [1,3,4];
nrnotransp := 2 ;
listnotransp := [3,4] ;
#  Initial conditions
#  vic: options to select initial guesses
# - 1: read from a file "initialconc.txt", which contains the columns 
#      z, conc(species1), conc(species2)....
# - 2: fixed concentration, defined in the array iniconc
# - 3: individual files for each species. listinput contains the species number (lenght: number of species), filenames are given in "file_in_names".inp, containing in column 1 the conc, in column 2 depth. 
# 
# example (only relevant input data is given for each option, assuming 17 species):
# NOTE that you have to provide something for all ncomp species in the arrays iniconc, listinput and file_in_names, EVEN IF YOU DON'T USE IT WITH THE OPTION YOU SELECTED
# vic:=1;
# vic:=2;iniconc:=[0.,2.,...]; 
# vic:=3;listinput:=[1,2,...,17];file_in_names:=[o2,no3,...,sp17];
vic := 2 ;
iniconc := [1.0e-6,1.0e-6,1.0e-7,1.0e-7,1.0e-6,1.0e-6] ;
listinput := [1,2,3,4,5,6] ;
file_in_names := [dummy1,dummy2,dummy3,dummy4,dummy5,dummy6] ;
#  Output
# noutput: number of species to be printed
# nroutput: number of rates to be printed
# listoutput: species number to print
# listroutput: rate number to print
# file_names: Respective file name for each of the species to print
# file_rnames: Respective file name for each of the rates to print
# time_iniout: First time (in years) for which a printout is requested
# time_intvout: time interval (in years) at which the printing is performed, starting from time_iniout
# 
# example:
# noutput:=4;
# listoutput:=[1,2,3,5];
# file_names:=[o2, so4, MnOx, hco3];
# time_iniout:=10.;
# time_intvout:= 100.;
noutput := 6 ;
nroutput := 5 ;
listoutput := [1,2,3,4,5,6] ;
listroutput := [1,2,3,4,5] ;
file_names := [ch2o,c_bioav,baco,bacs,o2,so4] ;
file_rnames := [xrate1,xrate2,xrate3,xrate4,xrate5] ;
time_iniout := 0.01 ;
time_intvout := 250. ;
#  Optimization
# there are several optimization options available. here you can specify what needs to be optimized and where the data is stored. to identify what kind of algorithm you want to use please select the appropriate options in drivervalues.f
# nopt_v: number of parameters to be optimized
# ntopt_v: number of time points where measurements are available
# nparam_opt: total number of parameters. can include the nparam but also the physical parameters.
# maxxmeas_v: maximum number of depth points at any given time measured (used to make array sizes)
# maxspmeas_v: maximum number of species measured at any given time (used to make array sizes)
# opt_name: names of the parameters. they have to match the names given above, so best you make a copy paste!
# idpar_v: identify the parameters to be optimized from the parameter list opt_name
# filemeas_name: name of the files containing the measured data
# if nopt_v is set to 0 then the rest of the input doesn't matter
# example:
# nopt_v := 2; # number of parameters to be optimized
# ntopt_v := 3; # number of timepoints with measurements
# nparam_opt := nparam; # total number of parameters, set equal to all parameters except physcial ones
# maxxmeas_v:=20; # maximum 20 points in a profile at any given time
# maxspmeas_v :=2; # maximum 2 species measured at one timepoint
# opt_name := bioname; # (note that this does not include the physical parameters! if you want them to be adapted you need to specify them explicitly)
# idpar_v:= [1,3]; # optimize parameters 1 and 3 in the above list
# filemeas_name:= [meas1.dat, meas2.dat, meas3.dat]; # filenames with measurements at timepoints
nopt_v := 0 ;
ntopt_v := 0 ;
nparam_opt := 0 ;
maxxmeas_v := 0;
maxspmeas_v := 0;
opt_name := [];
idpar_v := [];
filemeas_name := [];
#  Maple specific info
# dir_f: directory where the FORTRAN routines and Maple spread.m files are parsed
# format Mac: "Macinthosh HD:UU:...:code"
# format PC: "C:\\maple\\...\\code"
# WAS: dir_f := "C:\\Dokumente und Einstellungen\\centler\\Desktop\\Labor\\Simulations": 
# currentdir(dir_f):
# save "spread.m" ;
dir_f := "C:\\Program Files\\BRNSPackage" :
currentdir(dir_f) :
parse(sprintf("save %q,\"spread.m\";",anames()), statement) ;

"now execute processor - make sure the directories are set correctly";
# 
# ACG
