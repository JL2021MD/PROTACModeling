#!/bin/bash
#This script was ran on the cluster for each PMEMD folder.
cat >01.min.mdin <<EOF
Minmize all the hydrogens
&cntrl
 imin=1,           ! Minimize the initial structure
 ntmin=2,         ! Use steepest descent Ryota Added
 maxcyc=5000,    ! Maximum number of cycles for minimization
 igb=5,           !
 ntb=0,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 ntp=0,            ! No pressure scaling
 ntf=1,            ! Complete force evaluation
 ntwx= 1000,       ! Write to trajectory file every ntwx steps
 ntpr= 1000,       ! Print to mdout every ntpr steps
 ntwr= 1000,       ! Write a restart file every ntwr steps
 cut=999,         ! Nonbonded cutoff in Angstroms
 ntr=1,            ! Turn on restraints
 restraintmask=":1-513 & !@H=", ! atoms to be restrained (all in residue 1-513 but not H)
 restraint_wt=5.0, ! force constant for restraint
 ntxo=2,           ! Write coordinate file in NetCDF format (1 for ASCII, 2 for NetCDF)
 ioutfm=1,         ! Write trajectory file in NetCDF format (0 for ASCII)
/
EOF

cat >02.equil.mdin <<EOF
 MD simualation
&cntrl
 imin=0,           ! Perform MD
 nstlim=50000      ! Number of MD steps
 ntb=0,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 igb=5,            ! 1,2,5,8=implicit solvent
 ntc=1,            ! 2=SHAKE on for bonds with hydrogen (SHAKE restricts bond stretches) (3=on for all bonds, 1=off, see p352 manual)
 dt=0.001,         ! Timestep (ps)
 ntp=0,            ! 1=Isotropic pressure scaling (0=no scaling) (should be 1 or 2 when constant pressure periodic boundaries are used)
 barostat=2        ! 1=Berendsen barostat, 2=montecarlo thermostat, p339 amber manual, 2 is suppose to be faster when ntb=2
 taup=0.5          ! Pressure relaxtion time (ps)
 ntf=1,            ! 2=No force evaluation for bonds with hydrogen (1=complete force evaluation)
 ntt=3,            ! Langevin thermostat (0=constant total energy, p337 manual)
 gamma_ln=2.0      ! Collision Frequency for thermostat
 ig=-1,            ! Random seed for thermostat
 temp0=310         ! Simulation temperature (K)
 ntwx= 10000,        ! Write to trajectory file every ntwx steps
 ntpr= 1000,        ! Print to mdout every ntpr steps
 ntwr= 10000,        ! Write a restart file every ntwr steps
 cut=999,          ! Nonbonded cutoff in Angstroms
 ntr=1,            ! Turn on restraints
 restraintmask=":1-513 & !@H=", ! atoms to be restrained (all in residue 1-513 but not H)
 restraint_wt=5.0, ! force constant for restraint
 ntxo=2,           ! 1=Write coordinate file in ASCII format, 2=NetCDF
 ioutfm=1,         ! 0=Write trajectory file in ASCII format, 1=NetCDF
 nmropt=1,         !
/

&wt type='END' /
DISANG=output.rst
LISTOUT=POUT
EOF

cat>06.equil.mdin <<EOF
 MD simualation
&cntrl
 imin=0,           ! Perform MD
 nstlim=50000      ! Number of MD steps
 ntb=0,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 igb=5,            ! 1,2,5,8=implicit solvent
 ntc=1,            ! 2=SHAKE on for bonds with hydrogen (SHAKE restricts bond stretches) (3=on for all bonds, 1=off, see p352 manual)
 dt=0.001,         ! Timestep (ps)
 ntp=0,            ! 1=Isotropic pressure scaling (0=no scaling) (should be 1 or 2 when constant pressure periodic boundaries are used)
 barostat=1        ! 1=Berendsen barostat, 2=montecarlo thermostat, p339 amber manual, 2 is suppose to be faster when ntb=2
 taup=0.5          ! Pressure relaxtion time (ps)
 ntf=1,            ! 2=No force evaluation for bonds with hydrogen (1=complete force evaluation)
 ntt=3,            ! Langevin thermostat (0=constant total energy, p337 manual)
 gamma_ln=2.0      ! Collision Frequency for thermostat 
 ig=-1,            ! Random seed for thermostat
 temp0=310         ! Simulation temperature (K)
 ntwx= 10000,        ! Write to trajectory file every ntwx steps
 ntpr= 1000,        ! Print to mdout every ntpr steps
 ntwr= 10000,        ! Write a restart file every ntwr steps
 cut=999,          ! Nonbonded cutoff in Angstroms
 ntr=1,            ! Turn on restraints
 restraintmask=":128-512|:513@O1,O2,O3,O4,C22,C23,N5,C24,C25,C26,N4,C21,C20,C27,C28,C19,C18,C17,C16", ! atoms to be restrained
 restraint_wt=10, ! force constant for restraint
 ntxo=2,           ! 1=Write coordinate file in ASCII format, 2=NetCDF
 ioutfm=1,         ! 0=Write trajectory file in ASCII format, 1=NetCDF
 iwrap=0,          ! iwrap is turned on
 nmropt=1,         !
/
&wt type='END' /
DISANG=output.rst
LISTOUT=POUT
EOF
cat >07.equil.mdin <<EOF
 MD simualation
&cntrl
 imin=0,           ! Perform MD
 nstlim=50000      ! Number of MD steps
 ntb=0,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 igb=5,            ! 1,2,5,8=implicit solvent
 ntc=1,            ! No SHAKE on bonds between hydrogens
 dt=0.001,         ! Timestep (ps)
 ntp=0,            ! Isotropic pressure scaling
 barostat=1        ! Berendsen
 taup=0.5          ! Pressure relaxtion time (ps)
 ntf=1,            ! Complete force evaluation
 ntt=3,            ! Langevin thermostat
 gamma_ln=2.0      ! Collision Frequency for thermostat
 ig=-1,            ! Random seed for thermostat
 temp0=310         ! Simulation temperature (K)
 ntwx= 10000,       ! Write to trajectory file every ntwx steps
 ntpr= 1000,       ! Print to mdout every ntpr steps
 ntwr= 10000,       ! Write a restart file every ntwr steps
 cut=  999,        ! Nonbonded cutoff in Angstroms
 ntr=1,            ! Turn on restraints
 restraintmask=":128-512|:513@O1,O2,O3,O4,C22,C23,N5,C24,C25,C26,N4,C21,C20,C27,C28,C19,C18,C17,C16", ! atoms to be restrained
 restraint_wt=10, ! force constant for restraint
 ntxo=2,           ! 1=Write coordinate file in ASCII format, 2=NetCDF
 ioutfm=1,         ! 0=Write trajectory file in ASCII format, 1=NetCDF
 iwrap=0,          ! iwrap is turned off
 nmropt=1,         !
/
&wt type='END' /
DISANG=output.rst
LISTOUT=POUT

EOF

cat >08.equil.mdin <<EOF
MD simulations
&cntrl
 imin=0,           ! Perform MD
 nstlim=50000      ! Number of MD steps
 ntx=5,            ! Positions and velocities read formatted
 irest=1,          ! Restart calculation
 ntc=1,            ! No SHAKE on for bonds with hydrogen
 dt=0.001,         ! Timestep (ps)
 ntb=0,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 igb=5,            ! 1,2,5,8=implicit solvent
 ntp=0,            ! 1=Isotropic pressure scaling, 0=no scaling
 barostat=1        ! Berendsen
 taup=0.5          ! Pressure relaxtion time (ps)
 ntf=1,            ! Complete force evaluation
 ntt=3,            ! Langevin thermostat
 gamma_ln=2.0      ! Collision Frequency for thermostat
 ig=-1,            ! Random seed for thermostat
 temp0=310         ! Simulation temperature (K)
 ntwx= 10000,       ! Write to trajectory file every ntwx steps
 ntpr= 1000,       ! Print to mdout every ntpr steps
 ntwr= 10000,       ! Write a restart file every ntwr steps
 cut=  999,        ! Nonbonded cutoff in Angstroms
 ntr=1,            ! Turn on restraints
 restraintmask=":128-512|:513@O1,O2,O3,O4,C22,C23,N5,C24,C25,C26,N4,C21,C20,C27,C28,C19,C18,C17,C16", ! atoms to be restrained
 restraint_wt=10, ! force constant for restraint
 ntxo=2,           ! 1=Write coordinate file in ASCII format, 2=NetCDF
 ioutfm=1,         ! 0=Write trajectory file in ASCII format, 1=NetCDF
 iwrap=0,          ! iwrap is turned on
 nmropt=1,         !
/
&wt type='END' /
DISANG=output.rst
LISTOUT=POUT
EOF

cat >09.equil.mdin<<EOF
MD simulations
&cntrl
 imin=0,           ! Perform MD
 nstlim=50000      ! Number of MD steps
 ntx=5,            ! Positions and velocities read formatted
 irest=1,          ! Restart calculation
 ntc=1,            ! No SHAKE on for bonds with hydrogen
 dt=0.001,         ! Timestep (ps)
 ntb=0,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 igb=5,            ! 1,2,5,8=implicit solvent
 ntp=0,            ! 1=Isotropic pressure scaling,0=no scaling
 barostat=1        ! Berendsen
 taup=0.5          ! Pressure relaxtion time (ps)
 ntf=1,            ! Complete force evaluation
 ntt=3,            ! Langevin thermostat
 gamma_ln=2.0      ! Collision Frequency for thermostat
 ig=-1,            ! Random seed for thermostat
 temp0=310,         ! Simulation temperature (K)
 ntwx= 10000,       ! Write to trajectory file every ntwx steps
 ntpr= 1000,       ! Print to mdout every ntpr steps
 ntwr= 10000,       ! Write a restart file every ntwr steps
 cut=  999,        ! Nonbonded cutoff in Angstroms
 ntr=1,            ! Turn on restraints
 restraintmask=":128-512|:513@O1,O2,O3,O4,C22,C23,N5,C24,C25,C26,N4,C21,C20,C27,C28,C19,C18,C17,C16", ! atoms to be restrained
 restraint_wt=10, ! force constant for restraint
 ntxo=2,           ! 1=Write coordinate file in ASCII format
 ioutfm=1,         ! 0=Write trajectory file in ASCII format
 iwrap=0,          ! iwrap is turned on
 nmropt=1,         !
/
&wt type='END' /
DISANG=output.rst
LISTOUT=POUT
EOF

cat >10.prod.mdin <<EOF
 MD simulations
&cntrl
 imin=0,           ! Perform MD (1=energy minimization)
 nstlim=10000000   ! Number of MD steps
 ntx=5,            ! Both positions and velocities are read (1=only read positions)
 irest=1,          ! Restart ("continue" should be better word) calculation (0=start a new calculation from static)
 ntc=2,            ! SHAKE on for bonds with hydrogen (SHAKE restricts bond stretches) (3=on for all bonds, 1=off, see p352 manual)
 dt=0.002,         ! Timestep (ps)
 ntb=0,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 igb=5,            !
 ntp=0,            ! 1=Isotropic pressure scaling (0=no scaling) (should be 1 or 2 when constant pressure periodic boundaries are used)
 barostat=2        ! 1=Berendsen barostat, 2=montecarlo thermostat, p339 amber manual, 2 is suppose to be faster with ntb=2
 taup=0.5          ! Pressure relaxtion time (ps)
 ntf=2,            ! No force evaluation for bonds with hydrogen (1=complete force evaluation)
 ntt=3,            ! Langevin thermostat (0=constant total energy, p337 manual)
 gamma_ln=2.0      ! Collision Frequency for thermostat
 ig=-1,            ! Random number generator for thermostat (-1=based on current date and time)
 temp0=310,         ! Simulation temperature (K)
 ntwx= 50000,       ! Write to trajectory file every ntwx steps
 ntpr= 5000,       ! Print to mdout every ntpr steps
 ntwr= 50000,       ! Write a restart file every ntwr steps
 cut=  999,        ! Nonbonded cutoff in Angstroms
 ntr=1,            ! Turn on restraints 
 restraintmask=":128-512|:513@O1,O2,O3,O4,C22,C23,N5,C24,C25,C26,N4,C21,C20,C27,C28,C19,C18,C17,C16", ! atoms to be restrained
 restraint_wt=10, ! force constant for restraint
 ntxo=2,           ! Write coordinate file in NetCDF format (1=ASCII format)
 ioutfm=1,         ! Write trajectory file in NetCDF format (0=ASCII format)
 iwrap=0,          ! iwrap is turned on (0=off)
 nmropt=1,         !
/
&wt type='END' /
DISANG=output.rst
LISTOUT=POUT
EOF

source /home/amber20/amber.sh
echo "Started run on `date` "
do_cuda="pmemd.cuda" 

prmtop="complex.gas.prmtop"
coords="complex.gas" 


MDINPUTS=(01.min 02.equil 06.equil 07.equil 08.equil 09.equil 10.prod)

for input in ${MDINPUTS[@]}; do

 $do_cuda -O -i ${input}.mdin -o ${input}.mdout -p $prmtop -c ${coords}.rst7 -ref ${coords}.rst7 -x ${input}.trj -inf ${input}.info -r ${input}.rst7
 coords=$input
done 
rm *.mdin
echo "Finished run on `date` "
#######The following aligns each frame to the crystal pose##########
cat >cpptraj.align.in<<EOF
parm complex.gas.prmtop [1]
parm crystcomplex.gas.prmtop [2]


trajin complex.gas.rst7 parm [1]
trajin 02.equil.trj parm [1]
trajin 06.equil.trj parm [1]
trajin 07.equil.trj parm [1]
trajin 08.equil.trj parm [1]
trajin 09.equil.trj parm [1]
trajin 10.prod.trj parm [1]

reference crystcomplex.gas.rst7 parm [2]

rmsd rms1 reference :105@CB,CA:99@OD1:98@CG:41,42,44@CA:41@CG:40@CH2
trajout 2-10.pocket2.nc nobox
EOF
cpptraj cpptraj.align.in
#######The following generates RMSD values for different parts of the ternary complex to the crystal pose from the aligned frames##########
cat >RMSD.in<<EOF
parm crystcomplex.gas.prmtop [1]
parm complex.gas.prmtop [2]

trajin 2-10.pocket2.nc parm [2]

#read in reference                                                          
reference crystcomplex.gas.rst7 parm [1]
#compute rmsd and align CA to the BRD4 crystal structure (not starting structure)                        
rmsd BRD4toCrystal reference ":1-127@CA" out RMSDresults.csv nofit
rmsd CRBNtoCrystal reference ":129-512@CA" out RMSDresults.csv nofit
rmsd BRD4sidetoCrystal reference ":513@CL,C29,C41,C39,C,C2,C3,C5" out RMSDresults.csv nofit
rmsd CRBNsidetoCrystal reference ":513@C25,N5,C22,O3,O4,C19,O" out RMSDresults.csv nofit
EOF
cpptraj RMSD.in

#######The following calculates MMGBSA energy##########

#Define topology files

complex_prmtop="complex.gas.prmtop"
bothproteins_prmtop="2proteins.prmtop"
CRBN_prmtop="CRBN.prmtop"
BRD4_prmtop="BRD4.prmtop"
BRD4_PROTAC="BRD4+protac.prmtop"
CRBN_PROTAC="CRBN+protac.prmtop"
protac_prmtop="protac.prmtop"
trajectory="2-10.pocket2.nc"

#create mmgbsa input file
cat >mmgbsa.in<<EOF
mmgbsa CRBN-BRD4 analysis
&general                                                                                                                                   
  interval=1, netcdf=1, verbose=2
  keep_files=0, startframe=127, endframe=500
/                                                                                                                                          
&gb
  igb=5,                                          
  saltcon=0.0, surften=0.0072,
  surfoff=0.0, molsurf=0,
/
nmode 
  drms=0.002, maxcyc=10000,
  nminterval=1, nmendframe=2000,
  nmode_igb=1,                                                           
/                
                                               
EOF

mpirun -np 4 MMPBSA.py.MPI -O -i mmgbsa.in \
           -o  energy.results.dat \
           -eo energy.per-frame.dat \
           -cp ${complex_prmtop} \
           -rp ${bothproteins_prmtop} \
           -lp ${protac_prmtop} \
            -y ${trajectory}
sed -i '1i2 PROTEINS AS RECEPTOR
' energy.results.dat
sed -i '1i2 PROTEINS AS RECEPTOR
' energy.per-frame.dat

mpirun -np 4 MMPBSA.py.MPI -O -i mmgbsa.in \
           -o  BRD4_as_lig.dat \
           -eo BRD4_as_lig.per-frame.dat \
           -cp ${complex_prmtop} \
           -rp ${CRBN_PROTAC} \
           -lp ${BRD4_prmtop} \
            -y ${trajectory}
printf '\n\n\n\n\nBRD4 AS LIGAND \n' >> energy.results.dat
cat BRD4_as_lig.dat >> energy.results.dat

printf '\n\n\n\n\nBRD4 AS LIGAND \n' >> energy.per-frame.dat
cat BRD4_as_lig.per-frame.dat >> energy.per-frame.dat


mpirun -np 4 MMPBSA.py.MPI -O -i mmgbsa.in \
           -o  CRBN_as_lig.dat \
           -eo CRBN_as_lig.per-frame.dat \
           -cp ${complex_prmtop} \
           -rp ${BRD4_PROTAC} \
           -lp ${CRBN_prmtop} \
            -y ${trajectory}
printf '\n\n\n\n\nCRBN AS LIGAND \n' >> energy.results.dat
cat CRBN_as_lig.dat >> energy.results.dat

printf '\n\n\n\n\nCRBN AS LIGAND \n' >> energy.per-frame.dat
cat CRBN_as_lig.per-frame.dat >> energy.per-frame.dat

rm CRBN_as_lig.dat BRD4_as_lig.dat reference.frc CRBN_as_lig.per-frame.dat BRD4_as_lig.per-frame.dat mmgbsa.in
