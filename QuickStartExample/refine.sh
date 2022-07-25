#!/bin/bash
#run this script along with the inputs.in file in the folder to feed it the variables.
source /home/amber20/amber.sh
source inputs.in
ln -s ../Tleap/* .
cat > 01.min.mdin <<EOF
Minmize all the hydrogens
&cntrl
 imin=1,           ! 1=run minimzation, 0=run MD
 ntmin=1,         ! ncyc steps of steepest descent, followed by conjugate gradient min. 2=steepest descent only
 maxcyc=5000,    ! Total minimization steps
 ncyc=1000,      ! 1000 steps of steepest descent minimization
 igb=0,           ! 
 ntb=1,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 ntp=0,            ! No pressure scaling
 ntf=1,            ! Complete force evaluation
 ntwx= 1000,       ! Write to trajectory file every ntwx steps
 ntpr= 1000,       ! Print to mdout every ntpr steps
 ntwr= 1000,       ! Write a restart file every ntwr steps
 cut=8,         ! Nonbonded cutoff in Angstroms
 ntr=1,            ! Turn on restraints
 restraintmask="$min_eq_rst", ! atoms to be restrained (all heavy atoms but not H)
 restraint_wt=5.0, ! force constant for restraint
 ntxo=2,           ! Write coordinate file in NetCDF format (1 for ASCII, 2 for NetCDF)
 ioutfm=1,         ! Write trajectory file in NetCDF format (0 for ASCII)
/
EOF
cat > 02.equil.mdin<<EOF
 MD simualation
&cntrl
 imin=0,           ! 1=run minimzation, 0=run MD
 nstlim=50000      ! Number of MD steps
 ntb=2,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 igb=0,            ! 1,2,5,8=implicit solvent
 ntc=1,            ! 2=SHAKE on for bonds with hydrogen (SHAKE restricts bond stretches) (3=on for all bonds, 1=off, see p348 of 2020 manual)
 dt=0.001,         ! Timestep (ps)
 ntp=1,            ! 1=Isotropic pressure scaling (0=no scaling) (should be 1 or 2 when constant pressure periodic boundaries are used)
 barostat=2        ! 1=Berendsen barostat, 2=montecarlo thermostat, p347 of 2020 amber manual, 2 is suppose to be faster when ntb=2
 taup=0.5          ! Pressure relaxtion time (ps)
 ntf=1,            ! 2=No force evaluation for bonds with hydrogen (1=complete force evaluation)
 ntt=3,            ! Langevin thermostat (0=constant total energy, see p344 of 2020 manual)
 gamma_ln=2.0      ! Collision Frequency for thermostat
 ig=-1,            ! Random seed ftor thermostat
 temp0=310         ! Simulation temperature (K)
 ntwx= 10000,        ! Write to trajectory file every ntwx steps
 ntpr= 1000,        ! Print to mdout every ntpr steps
 ntwr= 10000,        ! Write a restart file every ntwr steps
 cut=8,          ! Nonbonded cutoff in Angstroms
 ntr=1,            ! Turn on restraints
 restraintmask="$min_eq_rst", ! atoms to be restrained
 restraint_wt=5.0, ! force constant for restraint
 ntxo=2,           ! 1=Write coordinate file in ASCII format, 2=NetCDF
 ioutfm=1,         ! 0=Write trajectory file in ASCII format, 1=NetCDF
 nmropt=0,         !
/
EOF
cat > 03.equil.mdin<<EOF
 MD simualation
&cntrl
 imin=0,           ! 1=run minimzation, 0=run MD
 nstlim=50000      ! Number of MD steps
 ntb=2,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 igb=0,            ! 1,2,5,8=implicit solvent
 ntc=1,            ! 2=SHAKE on for bonds with hydrogen (SHAKE restricts bond stretches) (3=on for all bonds, 1=off, see p348 of 2020 manual)
 dt=0.001,         ! Timestep (ps)
 ntp=1,            ! 1=Isotropic pressure scaling (0=no scaling) (should be 1 or 2 when constant pressure periodic boundaries are used)
 barostat=2        ! 1=Berendsen barostat, 2=montecarlo thermostat, p347 of 2020 amber manual, 2 is suppose to be faster when ntb=2
 taup=0.5          ! Pressure relaxtion time (ps)
 ntf=1,            ! 2=No force evaluation for bonds with hydrogen (1=complete force evaluation)
 ntt=3,            ! Langevin thermostat (0=constant total energy, see p344 of 2020 manual)
 gamma_ln=2.0      ! Collision Frequency for thermostat
 ig=-1,            ! Random seed for thermostat
 temp0=310         ! Simulation temperature (K)
 ntwx= 10000,        ! Write to trajectory file every ntwx steps
 ntpr= 1000,        ! Print to mdout every ntpr steps
 ntwr= 10000,        ! Write a restart file every ntwr steps
 cut=8,          ! Nonbonded cutoff in Angstroms
 ntr=1,            ! Turn on restraints
 restraintmask="$min_eq_rst", ! atoms to be restrained
 restraint_wt=1, ! force constant for restraint
 ntxo=2,           ! 1=Write coordinate file in ASCII format, 2=NetCDF
 ioutfm=1,         ! 0=Write trajectory file in ASCII format, 1=NetCDF
 nmropt=$nmropt,         !
/
&wt type='END' /
$more_rsts
EOF
cat > 04.equil.mdin<<EOF
 MD simualation
&cntrl
 imin=0,           ! 1=run minimzation, 0=run MD
 nstlim=50000      ! Number of MD steps
 ntb=2,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 igb=0,            ! 1,2,5,8=implicit solvent
 ntc=1,            ! 2=SHAKE on for bonds with hydrogen (SHAKE restricts bond stretches) (3=on for all bonds, 1=off, see p348 of 2020 manual)
 dt=0.001,         ! Timestep (ps)
 ntp=1,            ! 1=Isotropic pressure scaling (0=no scaling) (should be 1 or 2 when constant pressure periodic boundaries are used)
 barostat=2        ! 1=Berendsen barostat, 2=montecarlo thermostat, p347 of 2020 amber manual, 2 is suppose to be faster when ntb=2
 taup=0.5          ! Pressure relaxtion time (ps)
 ntf=1,            ! 2=No force evaluation for bonds with hydrogen (1=complete force evaluation)
 ntt=3,            ! Langevin thermostat (0=constant total energy, see p344 of 2020 manual)
 gamma_ln=2.0      ! Collision Frequency for thermostat
 ig=-1,            ! Random seed for thermostat
 temp0=310         ! Simulation temperature (K)
 ntwx= 10000,        ! Write to trajectory file every ntwx steps
 ntpr= 1000,        ! Print to mdout every ntpr steps
 ntwr= 10000,        ! Write a restart file every ntwr steps
 cut=8,          ! Nonbonded cutoff in Angstroms
 ntr=0,            ! Restraints OFF and next 2 lines are ineffective (1=on)
 restraintmask=":1-395 &!@H=", ! atoms to be restrained (all in residues 1-395 but not H)
 restraint_wt=10, ! force constant for restraint
 ntxo=2,           ! 1=Write coordinate file in ASCII format, 2=NetCDF
 ntxo=2,           ! 1=Write coordinate file in ASCII format, 2=NetCDF
 ioutfm=1,         ! 0=Write trajectory file in ASCII format, 1=NetCDF
 nmropt=$nmropt,         !
/
&wt type='END' /
$more_rsts

EOF
cat >10.prod.mdin<<EOF
 MD simulations
&cntrl
 imin=0,           ! Perform MD (1=energy minimization)
 nstlim=1000000    ! Number of MD steps
 ntx=5,            ! Both positions and velocities are read (1=only read positions)
 irest=1,          ! Continue calculation (0=start a new calculation from static)
 ntc=2,            ! 2=SHAKE on for bonds with hydrogen (SHAKE restricts bond stretches) (3=on for all bonds, 1=off, see p348 of 2020 manual)
 dt=0.004,         ! Timestep (ps)
 ntb=1,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 igb=0,            !
 ntp=0,            ! 1=Isotropic pressure scaling (0=no scaling) (should be 1 or 2 when constant pressure periodic boundaries are used)
 barostat=2        ! 1=Berendsen barostat, 2=montecarlo thermostat, p347 of 2020 amber manual, 2 is suppose to be faster when ntb=2
 taup=0.5          ! Pressure relaxtion time (ps)
 ntf=2,            ! No force evaluation for bonds with hydrogen (1=complete force evaluation)
 ntt=3,            ! Langevin thermostat (0=constant total energy, see p344 of 2020 manual)
 gamma_ln=0.01       ! Collision Frequency for thermostat
 ig=-1,            ! Random number generator for thermostat (-1=based on current date and time)
 temp0=310,         ! Simulation temperature (K)
 ntwx= 25000,       ! Write to trajectory file every ntwx steps
 ntpr= 5000,       ! Print to mdout every ntpr steps
 ntwr= 25000,       ! Write a restart file every ntwr steps
 cut=  8,        ! Nonbonded cutoff in Angstroms
 ntr=0,            ! Restraints OFF and next 2 lines are ineffective (1=on)
 restraintmask=":396>@8&!:WAT,Na+", ! atoms to be restrained (atoms 8A or more away from residue 396, but not Water or sodium)
 restraint_wt=10, ! force constant (Kcal/mol/A2) for restraint
 ntxo=2,           ! Write coordinate file in NetCDF format (1=ASCII format)
 ioutfm=1,         ! Write trajectory file in NetCDF format (0=ASCII format)
 iwrap=1,          ! iwrap is turned on (0=off)
 nmropt=$nmropt,         !
/
&wt type='END' /
$more_rsts

EOF
echo "Started Equil and MD on `date` "
do_cuda="pmemd.cuda" 

prmtop="complex.HMR.opc.prmtop"
coords="complex.opc" 

MDINPUTS=(01.min 02.equil 03.equil 04.equil 10.prod)

for input in ${MDINPUTS[@]}; do

 $do_cuda -O -i ${input}.mdin -o ${input}.mdout -p $prmtop -c ${coords}.rst7 -ref ${coords}.rst7 -x ${input}.trj -inf ${input}.info -r ${input}.rst7
 coords=$input
done 

echo "Finished MD on `date` "
rm *.mdin

cat > cpptraj.align.in <<EOF
parm complex.opc.prmtop [1]
parm complex.gas.prmtop [2]
trajin complex.opc.rst7
trajin 01.min.rst7
trajin 02.equil.rst7
trajin 03.equil.rst7
trajin 04.equil.rst7
trajin 10.prod.trj
strip :WAT:Na+:Cl-
autoimage

reference complex.gas.rst7 parm [2]

rms rms1 $aligned_forstrip reference
trajout 2-10.align.nc
trajout 2-10.pdb
EOF

cpptraj -i cpptraj.align.in
cat > RMSD.in <<EOF
parm complex.gas.prmtop
trajin 2-10.align.nc
reference complex.gas.rst7

rmsd $rmsd_unaligned_protein reference $rmsd_unaligned_residues out RMSDresults.csv nofit
rmsd $rmsd_aligned_protein reference $rmsd_aligned_residues out RMSDresults.csv nofit
rmsd $rmsd_unaligned_lig reference $rmsd_unaligned_ligatoms out RMSDresults.csv nofit
rmsd $rmsd_aligned_lig reference $rmsd_aligned_ligatoms out RMSDresults.csv nofit
rmsd linker reference $rmsd_linkeratoms out RMSDresults.csv nofit

EOF
cpptraj -i RMSD.in

##This saves your trajectory in PDB format (with pruned PDB files to reduce space), for the quickest visualization
awk '$0=substr($0,1,54)' 2-10.pdb > 2-10cut.pdb
tar -czf 2-10cut.pdb.tar.gz 2-10cut.pdb
rm 2-10cut.pdb 2-10.pdb

##Removing water trajectories to save space.
if [ -s RMSDresults.csv ]
then
rm 1_* docker* *trj *info 0*rst7 10*.rst7
fi

find -type l -delete
