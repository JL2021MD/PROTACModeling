#!/bin/bash
#This script was ran on the cluster for each PMEMD folder.
source /home/amber20/amber.sh
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
 restraintmask=":1-253 & !@H=", ! atoms to be restrained (all in residue 1-253 but not H)
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
 restraintmask=":1-253 & !@H=", ! atoms to be restrained (all in residue 1-253 but not H)
 restraint_wt=5.0, ! force constant for restraint
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
 ntr=0,            ! Turn on restraints
 restraintmask=":112-252 & !@H=|:253@SBR,CCF,NBM,OD1,CCN,CBD", ! atoms to be restrained (all in BRD4 but not H)
 restraint_wt=10, ! force constant for restraint
 ntxo=2,           ! 1=Write coordinate file in ASCII format, 2=NetCDF
 ioutfm=1,         ! 0=Write trajectory file in ASCII format, 1=NetCDF
 nmropt=0,         !
/

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
 ntr=0,            ! Turn off restraints
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
 nstlim=12500000   ! Number of MD steps
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
 gamma_ln=0.01      ! Collision Frequency for thermostat
 ig=-1,            ! Random number generator for thermostat (-1=based on current date and time)
 temp0=310,         ! Simulation temperature (K)
 ntwx= 25000,       ! Write to trajectory file every ntwx steps
 ntpr= 5000,       ! Print to mdout every ntpr steps
 ntwr= 25000,       ! Write a restart file every ntwr steps
 cut=  8,        ! Nonbonded cutoff in Angstroms
 ntr=0,            ! Turn off restraints 
 ntxo=2,           ! Write coordinate file in NetCDF format (1=ASCII format)
 ioutfm=1,         ! Write trajectory file in NetCDF format (0=ASCII format)
 iwrap=0,          ! iwrap is turned on (0=off)
 nmropt=0,         !
/

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
#######The following aligns each frame to the crystal pose##########
cat > cpptraj.align.in <<EOF
parm complex.opc.prmtop
parm crystcomplex.gas.prmtop [2]

trajin complex.opc.rst7
trajin 02.equil.trj
trajin 03.equil.trj
trajin 04.equil.trj
trajin 10.prod.trj
strip :WAT:Na+:Cl-

reference crystcomplex.gas.rst7 parm [2]

rms rms1 :7-107@CA reference
trajout 2-10.BRD4align.nc
EOF

cpptraj -i cpptraj.align.in
#######The following generates RMSD values for different parts of the ternary complex to the crystal pose from the aligned frames##########
cat > RMSD.in <<EOF
parm complex.gas.prmtop
parm crystcomplex.gas.prmtop [2]
trajin 2-10.BRD4align.nc
reference crystcomplex.gas.rst7 parm [2]

rmsd VHL reference ":112-252@CA" out RMSDresults.csv nofit
rmsd BRD4 reference ":1-111@CA" out RMSDresults.csv nofit
rmsd BRD4sidetoCrystal reference ":253@CL1,CBX,CAD,CAC,CAB,CCG,CCL,CBU" out RMSDresults.csv nofit
rmsd VHLsidetoCrystal reference ":253@CAA,CAV,CCF,CBC,C,OD1,N,CCQ,CBT" out RMSDresults.csv nofit
rmsd LinkertoCrystal reference ":253@OBO,CAY,CAZ,OBP,CBA,CBB,OBQ" out RMSDresults.csv nofit

EOF
cpptraj -i RMSD.in

#######The following calculates MMGBSA energy##########

complex_prmtop="complex.gas.prmtop"
bothproteins_prmtop="2proteins.prmtop"
VHL_prmtop="VHL.prmtop"
BRD4_prmtop="BRD4.prmtop"
BRD4_PROTAC="BRD4+protac.prmtop"
VHL_PROTAC="VHL+protac.prmtop"
protac_prmtop="protac.prmtop"
trajectory="2-10.BRD4align.nc"

#create mmgbsa input file
cat >mmgbsa.in<<EOF
mmgbsa VHL-BRD4 analysis
&general                                                                                                                                   
  interval=1, netcdf=1, verbose=2
  keep_files=0, startframe=417, endframe=1000
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

out="energy.results.dat"
out_frames="energy.per-frame.dat"

mpirun -np 4 MMPBSA.py.MPI -O -i mmgbsa.in \
           -o  $out \
           -eo $out_frames \
           -cp ${complex_prmtop} \
           -rp ${bothproteins_prmtop} \
           -lp ${protac_prmtop} \
            -y ${trajectory}
sed -i '1i2 PROTEINS AS RECEPTOR
' $out
sed -i '1i2 PROTEINS AS RECEPTOR
' $out_frames
rm reference.frc
mpirun -np 4 MMPBSA.py.MPI -O -i mmgbsa.in \
           -o  BRD4_as_lig.dat \
           -eo BRD4_as_lig.per-frame.dat \
           -cp ${complex_prmtop} \
           -rp ${VHL_PROTAC} \
           -lp ${BRD4_prmtop} \
            -y ${trajectory}
printf '\n\n\n\n\nBRD4 AS LIGAND \n' >> $out
cat BRD4_as_lig.dat >> $out

printf '\n\n\n\n\nBRD4 AS LIGAND \n' >> $out_frames
cat BRD4_as_lig.per-frame.dat >> $out_frames

rm reference.frc
mpirun -np 4 MMPBSA.py.MPI -O -i mmgbsa.in \
           -o  VHL_as_lig.dat \
           -eo VHL_as_lig.per-frame.dat \
           -cp ${complex_prmtop} \
           -rp ${BRD4_PROTAC} \
           -lp ${VHL_prmtop} \
            -y ${trajectory}
printf '\n\n\n\n\nVHL AS LIGAND \n' >> $out
cat VHL_as_lig.dat >> $out

printf '\n\n\n\n\nVHL AS LIGAND \n' >> $out_frames
cat VHL_as_lig.per-frame.dat >> $out_frames

rm VHL_as_lig.dat BRD4_as_lig.dat reference.frc VHL_as_lig.per-frame.dat BRD4_as_lig.per-frame.dat mmgbsa.in
