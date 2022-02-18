#!/bin/bash

ln -s ../Tleap/* .

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
 restraintmask=":1-267 & !@H=", ! atoms to be restrained (all in residue 1-253 but not H)
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
DISANG=VHLandLig10free1.rst
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
&wt type='END' /
DISANG=VHLandLig10free1.rst
LISTOUT=POUT

EOF
cat > 05.equil.dummy.mdin<<EOF
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
 temp0=300         ! Simulation temperature (K)
 ntwx= 10000,        ! Write to trajectory file every ntwx steps
 ntpr= 1000,        ! Print to mdout every ntpr steps
 ntwr= 10000,        ! Write a restart file every ntwr steps
 cut=999,          ! Nonbonded cutoff in Angstroms
 ntr=1,            ! Turn on restraints
 restraintmask=":1-163,176-512 & !@H=", ! atoms to be restrained (all in residue 1-513 but not H)
 restraint_wt=0.1, ! force constant for restraint
 ntxo=2,           ! 1=Write coordinate file in ASCII format, 2=NetCDF
 ioutfm=1,         ! 0=Write trajectory file in ASCII format, 1=NetCDF
 nmropt=1,         !
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
&wt type='TEMP0', istep1=0,istep2=25000000,
   value1=310.0, value2=410.0
/
&wt type='END'
/

EOF

cat >11.prod.mdin<<EOF
 MD simulations
&cntrl
 imin=0,           ! Perform MD (1=energy minimization)
 nstlim=45000000   ! Number of MD steps
 ntx=5,            ! Both positions and velocities are read (1=only read positions)
 irest=1,          ! Restart ("continue" should be better word) calculation (0=start a new calculation from static)
 ntc=2,            ! SHAKE on for bonds with hydrogen (SHAKE restricts bond stretches) (3=on for all bonds, 1=off, see p352 manual)
 dt=0.002,         ! Timestep (ps)
 ntb=2,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 igb=0,            !
 ntp=1,            ! 1=Isotropic pressure scaling (0=no scaling) (should be 1 or 2 when constant pressure periodic boundaries are used)
 barostat=2        ! 1=Berendsen barostat, 2=montecarlo thermostat, p339 amber manual, 2 is suppose to be faster with ntb=2
 taup=0.5          ! Pressure relaxtion time (ps)
 ntf=2,            ! No force evaluation for bonds with hydrogen (1=complete force evaluation)
 ntt=3,            ! Langevin thermostat (0=constant total energy, p337 manual)
 gamma_ln=0.01     ! Collision Frequency for thermostat
 ig=-1,            ! Random number generator for thermostat (-1=based on current date and time)
 temp0=310,         ! Simulation temperature (K)
 ntwx= 50000,       ! Write to trajectory file every ntwx steps
 ntpr= 5000,       ! Print to mdout every ntpr steps
 ntwr= 50000,       ! Write a restart file every ntwr steps
 cut=  8,        ! Nonbonded cutoff in Angstroms
 ntr=0,            ! Turn on restraints 
 restraintmask=":1-111 & !@H=|:253@CCE,CCL,CAB,CCJ", ! atoms to be restrained (all in residue 1-513 but not H)
 restraint_wt=10, ! force constant for restraint
 ntxo=2,           ! Write coordinate file in NetCDF format (1=ASCII format)
 ioutfm=1,         ! Write trajectory file in NetCDF format (0=ASCII format)
 iwrap=0,          ! iwrap is turned on (0=off)
 nmropt=1,         !
/
&wt type='TEMP0', istep1=0,istep2=45000000,
   value1=310.0, value2=400.0
/
&wt type='END'
/

EOF

source /home/amber20/amber.sh
#if false; then
echo "Started Equil and MD on `date` "
do_cuda="pmemd.cuda" 

prmtop="complex.HMR.opc.prmtop"
coords="04.min" 


MDINPUTS=(02.equil 03.equil 04.equil 10.prod)
#MDINPUTS=(10.prod)

for input in ${MDINPUTS[@]}; do

# $do_cuda -O -i ${input}.mdin -o ${input}.mdout -p $prmtop -c ${coords}.rst7 -ref ${coords}.rst7 -x ${input}.trj -inf ${input}.info -r ${input}.rst7
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
rms rms1 :5-116@CA reference
trajout 2-10.BDalign.nc
trajout 2-10.pdb
EOF

cpptraj -i cpptraj.align.in
cat > RMSD.in <<EOF
parm complex.gas.prmtop
trajin 2-10.BDalign.nc
reference complex.gas.pdb

rmsd VHL reference ":118-209,251-266@CA" out RMSDresults.csv nofit
rmsd BRD4 reference ":5-116@CA" out RMSDresults.csv nofit
#rmsd BRD4sidetoCrystal reference ":267@O6,C43,C39,C40,N6,N5" out RMSDresults.csv nofit
#rmsd VHLsidetoCrystal reference ":267@C3,C4,C8,C6,C10,O2,C15,C12,N4" out RMSDresults.csv nofit
#rmsd linker reference ":267@C38,O39,C41,C44,C61,C64,N41" out RMSDresults.csv nofit

EOF
#cpptraj -i RMSD.in
cat > RMSDstart.in <<EOF
parm complex.gas.prmtop
trajin 2-10.BDalign.nc
reference complex.gas.pdb

rms rms1 :5-116@CA reference

rmsd VHL reference ":118-209,251-266@CA" out RMSDstart.csv nofit
rmsd BRD4 reference ":5-116@CA" out RMSDstart.csv nofit
rmsd BRD4sidetoCrystal reference ":267@O6,C43,C39,C40,N6,N5" out RMSDstart.csv nofit
rmsd VHLsidetoCrystal reference ":267@C3,C4,C8,C6,C10,O2,C15,C12,N4" out RMSDstart.csv nofit
EOF
cpptraj -i RMSDstart.in
exit
#fi 
cat >mmgbsa.in<<EOF
mmgbsa VHL-BRD4 analysis
&general                                                                                                                                   
  interval=1, netcdf=1, verbose=2
  keep_files=0, startframe=6, endframe=1000
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

complex_prmtop="complex.gas.prmtop"
out=(complex.en BDasLig.en VHLasLig.en)
out_frames=(complex.per-frame.dat BDasLig.per-frame.dat VHLasLig.per-frame.dat)
receptor_prmtop=(2proteins.prmtop VHL+protac.prmtop BD+protac.prmtop)
ligand_prmtop=(protac.prmtop BD.prmtop VHL.prmtop)
trajectory="2-10.BDalign.nc"
for ((i=0; i < 3 ; i++ )); do

mpirun -np 4 MMPBSA.py.MPI -O -i mmgbsa.in \
           -o  ${out[i]} \
           -eo ${out_frames[i]} \
           -cp ${complex_prmtop} \
           -rp ${receptor_prmtop[i]} \
           -lp ${ligand_prmtop[i]} \
            -y ${trajectory}

printf "\n${out[i]}\n" >> energy.results.dat
printf "\n\n\n\n\n${out_frames[i]}\n" >> energy.perframe.dat
cat ${out[i]} >> energy.results.dat
cat ${out_frames[i]} >> energy.perframe.dat

done
rm *en reference.frc *per-frame.dat mmgbsa.in temp.temp
find -type l -delete
exit

