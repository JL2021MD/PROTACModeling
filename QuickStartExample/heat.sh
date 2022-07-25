#!/bin/bash
#run this script along with the inputs.in file in the folder to feed it the variables.
ln -s ../Tleap/* .
source /home/amber20/amber.sh
source inputs.in
cat > 01.min.mdin <<EOF
Minmize all the hydrogens
&cntrl
 imin=1,           ! Minimize the initial structure
  ntmin=1,         ! ncyc steps of steepest descent, followed by conjugate gradient min. 2=steepest descent only
 maxcyc=5000,    ! Maximum number of cycles for minimization
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
 restraintmask="$min_eq_rst", ! atoms to be restrained (all in residue 1-253 but not H)
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
 ntc=1,            ! 2=SHAKE on for bonds with hydrogen (SHAKE restricts bond stretches) (3=on for all bonds, 1=off, see p352 manual)
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
 restraintmask="$min_eq_rst", ! atoms to be restrained (all in residue 1-253 but not H)
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
 ntc=1,            ! 2=SHAKE on for bonds with hydrogen (SHAKE restricts bond stretches) (3=on for all bonds, 1=off, see p352 manual)
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
 restraintmask="$min_eq_rst", ! atoms to be restrained (all in BRD4 but not H)
 restraint_wt=1, ! force constant for restraint
 ntxo=2,           ! 1=Write coordinate file in ASCII format, 2=NetCDF
 ioutfm=1,         ! 0=Write trajectory file in ASCII format, 1=NetCDF
 nmropt=$nmropt,         !
/
&wt type='END'
/
$more_rsts
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
 ntr=0,            ! Restraints OFF and next 2 lines are ineffective (1=on)
 restraintmask=":1-111 & !@H=|:253@CCE,CCL,CAB,CCJ", ! atoms to be restrained (all in BRD4 but not H)
 restraint_wt=10, ! force constant for restraint
 ntxo=2,           ! 1=Write coordinate file in ASCII format, 2=NetCDF
 ntxo=2,           ! 1=Write coordinate file in ASCII format, 2=NetCDF
 ioutfm=1,         ! 0=Write trajectory file in ASCII format, 1=NetCDF
 nmropt=$nmropt,         !
/
&wt type='END'
/
$more_rsts

EOF

cat >10.prod.mdin<<EOF
 MD simulations
&cntrl
 imin=0,           ! Perform MD (1=energy minimization)
 nstlim=12500000  ! Number of MD steps
 ntx=5,            ! Both positions and velocities are read (1=only read positions)
 irest=1,          ! Restart ("continue" should be better word) calculation (0=start a new calculation from static)
 ntc=2,            ! SHAKE on for bonds with hydrogen (SHAKE restricts bond stretches) (3=on for all bonds, 1=off, see p352 manual)
 dt=0.002,         ! Timestep (ps)
 ntb=2,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 igb=0,            !
 ntp=1,            ! 1=Isotropic pressure scaling (0=no scaling) (should be 1 or 2 when constant pressure periodic boundaries are used)
 barostat=2        ! 1=Berendsen barostat, 2=montecarlo thermostat, p347 of 2020 amber manual, 2 is suppose to be faster when ntb=2
 taup=0.5          ! Pressure relaxtion time (ps)
 ntf=2,            ! No force evaluation for bonds with hydrogen (1=complete force evaluation)
 ntt=3,            ! Langevin thermostat (0=constant total energy, p337 manual)
 gamma_ln=0.01     ! Collision Frequency for thermostat
 ig=-1,            ! Random number generator for thermostat (-1=based on current date and time)
 temp0=310,         ! Simulation temperature (K)
 ntwx= 25000,       ! Write to trajectory file every ntwx steps
 ntpr= 5000,       ! Print to mdout every ntpr steps
 ntwr= 25000,       ! Write a restart file every ntwr steps
 cut=  8,        ! Nonbonded cutoff in Angstroms
 ntr=0,            ! Restraints OFF and next 2 lines are ineffective (1=on)
 restraint_wt=10, ! force constant for restraint
 ntxo=2,           ! Write coordinate file in NetCDF format (1=ASCII format)
 ioutfm=1,         ! Write trajectory file in NetCDF format (0=ASCII format)
 iwrap=1,          ! iwrap is turned on (0=off)
 nmropt=$nmropt,   !
/
&wt type='TEMP0', istep1=0,istep2=12500000,
   value1=310.0, value2=410.0
/
&wt type='END'
/
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
trajin 02.equil.trj
trajin 03.equil.trj
trajin 04.equil.trj
trajin 10.prod.trj

strip :WAT:Na+:Cl-
autoimage

reference complex.gas.rst7 parm [2]

rms rms1 $aligned_forstrip reference
trajout 2-10.align.nc
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

##Removing water trajectories to save space.
if [ -s RMSDresults.csv ]
then
rm 1_* docker* *trj *info 0*rst7 10*.rst7
fi

find -type l -delete





## This measures the occupancy and temperature scores for each single trajectory.
cat>HAPOD.py<<EOF
import sys

occuThres=3.5
disThres=4
rejoinThres=3
RMSDcol=3
#First column is 0
rawFrames={}
finalFrames={}
file=open(sys.argv[1],"r")
output=open('HAPOD.out', "w")
lines=file.readlines()

for line in lines[1:]:
    frame=int(line.split()[0])
    RMSD=float(line.split()[RMSDcol])
    rawFrames[frame]=RMSD

MoveAvgFrames={}
a=int(list(rawFrames)[-1])
print("A total of ",a, "frames read.", file=output)
print("A total of ",a, "frames read.")
for i in range(10,a+1):
    newFrame=0
    for k in range(0,10):
        one=(rawFrames[i-k])
        newFrame=one+newFrame

    outFrame=round(newFrame/10,4)
    finalFrames[i-9]=outFrame

##Because the first 16 frames are part of the equilibrum only not yet started any heating, but the moving average of 10 still counts them.
rawOccupany=-7
for i in range(1,a-8):
    if finalFrames[i]<occuThres:
        rawOccupany=rawOccupany+1
occupancy=rawOccupany/5
print("Occupancy Score: ", occupancy, file=output)
print("Occupancy Score: ", occupancy)

dissociation={}
counter=0
dis=0
k=0
print("Dissociation frames:", file=output)
print("Dissociation frames:")
for i in range(1,a-8):
    if finalFrames[i]>disThres and dis==0:
        k=i+k
        print(i, file=output)
        print(i)
        dis=1
        counter=counter+1
    elif finalFrames[i]<rejoinThres and dis==1:
        dis=0
    if i==a-9 and finalFrames[a-9] <disThres and dis==0:
        k=k+562
        print(562, file=output)
        print(562)
        counter=counter+1

dissocation=round((k/counter-12)/5+310,1)
print("Temperature Score: ",dissocation, file=output)
print("Temperature Score: ",dissocation)

EOF
python3 HAPOD.py RMSDresults.csv




