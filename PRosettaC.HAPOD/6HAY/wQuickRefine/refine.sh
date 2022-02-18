#!/bin/bash

ln -s ../Tleap/* .

cat > 01.min.mdin <<EOF
Minmize all the hydrogens
&cntrl
 imin=1,           ! Minimize the initial structure
 ntmin=2,         ! Use steepest descent Ryota Added
 maxcyc=5000,    ! Maximum number of cycles for minimization
 igb=0,           !
 ntb=1,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 ntp=0,            ! No pressure scaling
 ntf=1,            ! Complete force evaluation
 ntwx= 1000,       ! Write to trajectory file every ntwx steps
 ntpr= 1000,       ! Print to mdout every ntpr steps
 ntwr= 1000,       ! Write a restart file every ntwr steps
 cut=8,         ! Nonbonded cutoff in Angstroms
 ntr=1,            ! Turn on restraints
 restraintmask=":1-273 & !@H=", ! atoms to be restrained (all in residue 1-253 but not H)
 restraint_wt=5.0, ! force constant for restraint
 ntxo=2,           ! Write coordinate file in NetCDF format (1 for ASCII, 2 for NetCDF)
 ioutfm=1,         ! Write trajectory file in NetCDF format (0 for ASCII)
/
EOF

cat > 02.equil.mdin<<EOF
 MD simualation
&cntrl
 imin=0,           ! Perform MD
 nstlim=50000      ! Number of MD steps
 ntb=2,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 igb=0,            ! 1,2,5,8=implicit solvent
 ntc=1,            ! 2=SHAKE on for bonds with hydrogen (SHAKE restricts bond stretches) (3=on for all bonds, 1=off, see p352 manual)
 dt=0.001,         ! Timestep (ps)
 ntp=1,            ! 1=Isotropic pressure scaling (0=no scaling) (should be 1 or 2 when constant pressure periodic boundaries are used)
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
 cut=8,          ! Nonbonded cutoff in Angstroms
 ntr=1,            ! Turn on restraints
 restraintmask=":1-273 & !@H=", ! atoms to be restrained (all in residue 1-253 but not H)
 restraint_wt=2.0, ! force constant for restraint
 ntxo=2,           ! 1=Write coordinate file in ASCII format, 2=NetCDF
 ioutfm=1,         ! 0=Write trajectory file in ASCII format, 1=NetCDF
 nmropt=0,         !
/

EOF
cat > 03.equil.mdin<<EOF
 MD simualation
&cntrl
 imin=0,           ! Perform MD
 nstlim=50000      ! Number of MD steps
 ntb=2,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 igb=0,            ! 1,2,5,8=implicit solvent
 ntc=1,            ! 2=SHAKE on for bonds with hydrogen (SHAKE restricts bond stretches) (3=on for all bonds, 1=off, see p352 manual)
 dt=0.001,         ! Timestep (ps)
 ntp=1,            ! 1=Isotropic pressure scaling (0=no scaling) (should be 1 or 2 when constant pressure periodic boundaries are used)
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
 cut=8,          ! Nonbonded cutoff in Angstroms
 ntxo=2,           ! 1=Write coordinate file in ASCII format, 2=NetCDF
 ioutfm=1,         ! 0=Write trajectory file in ASCII format, 1=NetCDF
 nmropt=0,         !
/
&wt type='END' /
DISANG=ligs.rst
LISTOUT=POUT
EOF
cat > 04.equil.mdin<<EOF
 MD simualation
&cntrl
 imin=0,           ! Perform MD
 nstlim=50000      ! Number of MD steps
 ntb=2,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 igb=0,            ! 1,2,5,8=implicit solvent
 ntc=1,            ! 2=SHAKE on for bonds with hydrogen (SHAKE restricts bond stretches) (3=on for all bonds, 1=off, see p352 manual)
 dt=0.001,         ! Timestep (ps)
 ntp=1,            ! 1=Isotropic pressure scaling (0=no scaling) (should be 1 or 2 when constant pressure periodic boundaries are used)
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
 cut=8,          ! Nonbonded cutoff in Angstroms
 ntr=0,            ! Turn on restraints
 restraintmask=":1-111 & !@H=|:253@CCE,CCL,CAB,CCJ", ! atoms to be restrained (all in BRD4 but not H)
 restraint_wt=10, ! force constant for restraint
 ntxo=2,           ! 1=Write coordinate file in ASCII format, 2=NetCDF
 ntxo=2,           ! 1=Write coordinate file in ASCII format, 2=NetCDF
 ioutfm=1,         ! 0=Write trajectory file in ASCII format, 1=NetCDF
 nmropt=0,         !
/

EOF

cat >10.prod.mdin<<EOF
 MD simulations
&cntrl
 imin=0,           ! Perform MD (1=energy minimization)
 nstlim=1000000   ! Number of MD steps
 ntx=5,            ! Both positions and velocities are read (1=only read positions)
 irest=1,          ! Restart ("continue" should be better word) calculation (0=start a new calculation from static)
 ntc=2,            ! SHAKE on for bonds with hydrogen (SHAKE restricts bond stretches) (3=on for all bonds, 1=off, see p352 manual)
 dt=0.004,         ! Timestep (ps)
 ntb=1,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 igb=0,            !
 ntp=0,            ! 1=Isotropic pressure scaling (0=no scaling) (should be 1 or 2 when constant pressure periodic boundaries are used)
 barostat=2        ! 1=Berendsen barostat, 2=montecarlo thermostat, p339 amber manual, 2 is suppose to be faster with ntb=2
 taup=0.5          ! Pressure relaxtion time (ps)
 ntf=2,            ! No force evaluation for bonds with hydrogen (1=complete force evaluation)
 ntt=3,            ! Langevin thermostat (0=constant total energy, p337 manual)
 gamma_ln=0.01       ! Collision Frequency for thermostat
 ig=-1,            ! Random number generator for thermostat (-1=based on current date and time)
 temp0=310,         ! Simulation temperature (K)
 ntwx= 25000,       ! Write to trajectory file every ntwx steps
 ntpr= 5000,       ! Print to mdout every ntpr steps
 ntwr= 25000,       ! Write a restart file every ntwr steps
 cut=  8,        ! Nonbonded cutoff in Angstroms
 ntr=0,            ! Turn on restraints 
 restraintmask=":1-111 & !@H=|:253@CCE,CCL,CAB,CCJ", ! atoms to be restrained (all in residue 1-513 but not H)
 restraint_wt=10, ! force constant for restraint
 ntxo=2,           ! Write coordinate file in NetCDF format (1=ASCII format)
 ioutfm=1,         ! Write trajectory file in NetCDF format (0=ASCII format)
 iwrap=1,          ! iwrap is turned on (0=off)
 nmropt=0,         !
/

EOF


source /home/amber20/amber.sh

echo "Started Equil and MD on `date` "
do_cuda="pmemd.cuda" 

prmtop="complex.HMR.opc.prmtop"
coords="04.min" 


MDINPUTS=(02.equil 03.equil 04.equil 10.prod)


for input in ${MDINPUTS[@]}; do

  $do_cuda -O -i ${input}.mdin -o ${input}.mdout -p $prmtop -c ${coords}.rst7 -ref ${coords}.rst7 -x ${input}.trj -inf ${input}.info -r ${input}.rst7
 coords=$input
done 

echo "Finished MD on `date` "
rm *.mdin

cat > cpptraj.align.in <<EOF
parm complex.HMR.opc.prmtop [1]
parm complex.gas.pdb [2]
trajin complex.opc.rst7 
trajin 04.min.rst7
trajin 02.equil.rst7
trajin 03.equil.rst7
trajin 04.equil.rst7
trajin 10.prod.trj

strip :WAT:Na+:Cl-
autoimage

reference complex.gas.pdb parm [2]
#rms rms1 :27,37,51,91@CA reference
rms rms1 :7-120@CA reference
trajout 2-10.BRD4align.nc
trajout 2-10.pdb
EOF

cpptraj -i cpptraj.align.in
cat > RMSD.in <<EOF
parm complex.gas.prmtop [1]
trajin 2-10.BRD4align.nc

reference complex.gas.rst7 parm [1]

rmsd VHL reference ":121-213,256-272@CA" out RMSDresults.csv nofit
rmsd VHLsidetoCrystal reference ":273@C3,C4,C6,C8,C10,C12,O2,N4" out RMSDresults.csv nofit
rmsd linker reference ":273@C30,C29,O6,C38" out RMSDresults.csv nofit
rmsd BRD4 reference ":1-120@CA" out RMSDresults.csv nofit
rmsd BRD4sidetoCrystal reference ":273@O8,C42,C38,C39,N6,N5" out RMSDresults.csv nofit

EOF
cpptraj -i RMSD.in

find -type l -delete
