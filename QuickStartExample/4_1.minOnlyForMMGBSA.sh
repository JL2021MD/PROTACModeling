#!/bin/bash
for i in {1..4}; do
cd $i/
mkdir formin
cd formin

#The total run time with a mid-level GPU is expected to be about 10 minutes for 4 minimizations + 4 single-point MMGBSA calculations, which is run in serial. Pose 5 will need a CPU version of pmemd to start the minimization (see line 156).
ln -s ../Tleap/* .
ln -s ../../inputs.in .
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
 restraint_wt=1.0, ! force constant for restraint
 ntxo=2,           ! Write coordinate file in NetCDF format (1 for ASCII, 2 for NetCDF)
 ioutfm=1,         ! Write trajectory file in NetCDF format (0 for ASCII)
/
EOF
cat > 02.min.mdin<<EOF
 MD simualation
&cntrl
 imin=1,           ! 1=run minimzation, 0=run MD
 maxcyc=10000,    ! Maximum number of cycles for minimization
 ncyc=1000,      ! 1000 steps of steepest descent minimization
 ntb=1,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 igb=0,            ! 1,2,5,8=implicit solvent
 ntc=1,            ! 2=SHAKE on for bonds with hydrogen (SHAKE restricts bond stretches) (3=on for all bonds, 1=off, see p352 manual)
 dt=0.001,         ! Timestep (ps)
 ntp=0,            ! 1=Isotropic pressure scaling (0=no scaling) (should be 1 or 2 when constant pressure periodic boundaries are used)
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
 ntr=0,            ! Turn on restraints
 restraintmask="$min_eq_rst", ! atoms to be restrained (all in residue 1-253 but not H)
 restraint_wt=5.0, ! force constant for restraint
 ntxo=2,           ! 1=Write coordinate file in ASCII format, 2=NetCDF
 ioutfm=1,         ! 0=Write trajectory file in ASCII format, 1=NetCDF
 nmropt=0,         !
/
EOF


echo "Started Equil and MD on `date` "
do_cuda="pmemd.cuda" 
do_cpu="mpirun -np 16 pmemd.MPI"
prmtop="complex.HMR.opc.prmtop"
coords="complex.opc" 

MDINPUTS=(01.min )

#just in case the CPU code is needed, change the do_cuda to do_cpu
for input in ${MDINPUTS[@]}; do

 $do_cuda -O -i ${input}.mdin -o ${input}.mdout -p $prmtop -c ${coords}.rst7 -ref ${coords}.rst7 -x ${input}.trj -inf ${input}.info -r ${input}.rst7
 coords=$input
done 

input="02.min"
$do_cuda -O -i ${input}.mdin -o ${input}.mdout -p $prmtop -c ${coords}.rst7 -ref ${coords}.rst7 -x ${input}.trj -inf ${input}.info -r ${input}.rst7
 coords=$input

echo "Finished MD on `date` "
rm *.mdin *.trj *.info 

cat > cpptraj.align.in <<EOF
parm complex.opc.prmtop [1]
parm complex.gas.prmtop [2]
trajin complex.opc.rst7
trajin 01.min.rst7
trajin 02.min.rst7

strip :WAT:Na+:Cl-
autoimage

reference complex.gas.rst7 parm [2]

rms rms1 $aligned_forstrip reference
trajout 0-2.align.nc
trajout 0-2.pdb
EOF

cpptraj -i cpptraj.align.in


cat >mmgbsa.in<<EOF
mmgbsa CRBN-BRD4 analysis
&general
  interval=1, netcdf=1, verbose=2
  keep_files=0, startframe=3, endframe=3
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
out=(complex.en BRD4asLig.en CRBNasLig.en)
out_frames=(complex.per-frame.dat BRD4asLig.per-frame.dat CRBNasLig.per-frame.dat)
receptor_prmtop=(2proteins.prmtop CRBN+protac.prmtop BRD4+protac.prmtop)
ligand_prmtop=(protac.prmtop BRD4.prmtop CRBN.prmtop)
trajectory="0-2.align.nc"
for ((i=0; i < 3 ; i++ )); do

mpirun -np 1 MMPBSA.py.MPI -O -i mmgbsa.in \
           -o  ${out[i]} \
           -eo ${out_frames[i]} \
           -cp ${complex_prmtop} \
           -rp ${receptor_prmtop[i]} \
           -lp ${ligand_prmtop[i]} \
            -y ${trajectory}

printf "\n${out[i]}\n" >> energy.dat
printf "\n\n\n\n\n${out_frames[i]}\n" >> energy.perframe.dat
cat ${out[i]} >> energy.dat
cat ${out_frames[i]} >> energy.perframe.dat

done
rm *en reference.frc *per-frame.dat mmgbsa.in
find -type l -delete


cd ../../
done

## For pose 5 which needs a CPU pmemd start, we'll save a bit of time to link over the 01.min.rst7

cd 5/
mkdir formin
cd formin

ln -s ../Tleap/* .
ln -s ../../inputs.in .
ln -s ../../Other/MinOnly.Pose5.Refine/01.min.rst7 .
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
 restraint_wt=1.0, ! force constant for restraint
 ntxo=2,           ! Write coordinate file in NetCDF format (1 for ASCII, 2 for NetCDF)
 ioutfm=1,         ! Write trajectory file in NetCDF format (0 for ASCII)
/
EOF
cat > 02.min.mdin<<EOF
 MD simualation
&cntrl
 imin=1,           ! 1=run minimzation, 0=run MD
 maxcyc=10000,    ! Maximum number of cycles for minimization
 ncyc=1000,      ! 1000 steps of steepest descent minimization
 ntb=1,            ! 2=Constant Pressure, default when ntp>0. Unit cell size changes to keep P constant) (0=no periodicity applied and PME is off, default when igb>0)(1=constant volume, default when igb and ntp both 0)
 igb=0,            ! 1,2,5,8=implicit solvent
 ntc=1,            ! 2=SHAKE on for bonds with hydrogen (SHAKE restricts bond stretches) (3=on for all bonds, 1=off, see p352 manual)
 dt=0.001,         ! Timestep (ps)
 ntp=0,            ! 1=Isotropic pressure scaling (0=no scaling) (should be 1 or 2 when constant pressure periodic boundaries are used)
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
 ntr=0,            ! Turn on restraints
 restraintmask="$min_eq_rst", ! atoms to be restrained (all in residue 1-253 but not H)
 restraint_wt=5.0, ! force constant for restraint
 ntxo=2,           ! 1=Write coordinate file in ASCII format, 2=NetCDF
 ioutfm=1,         ! 0=Write trajectory file in ASCII format, 1=NetCDF
 nmropt=0,         !
/
EOF


echo "Started Equil and MD on `date` "
do_cuda="pmemd.cuda" 
do_cpu="mpirun -np 16 pmemd.MPI"
prmtop="complex.HMR.opc.prmtop"
coords="complex.opc" 

input="02.min"
$do_cuda -O -i ${input}.mdin -o ${input}.mdout -p $prmtop -c ${coords}.rst7 -ref ${coords}.rst7 -x ${input}.trj -inf ${input}.info -r ${input}.rst7
 coords=$input

echo "Finished MD on `date` "
rm *.mdin *.trj *.info 

cat > cpptraj.align.in <<EOF
parm complex.opc.prmtop [1]
parm complex.gas.prmtop [2]
trajin complex.opc.rst7
trajin 01.min.rst7
trajin 02.min.rst7

strip :WAT:Na+:Cl-
autoimage

reference complex.gas.rst7 parm [2]

rms rms1 $aligned_forstrip reference
trajout 0-2.align.nc
trajout 0-2.pdb
EOF

cpptraj -i cpptraj.align.in


cat >mmgbsa.in<<EOF
mmgbsa CRBN-BRD4 analysis
&general
  interval=1, netcdf=1, verbose=2
  keep_files=0, startframe=3, endframe=3
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
out=(complex.en BRD4asLig.en CRBNasLig.en)
out_frames=(complex.per-frame.dat BRD4asLig.per-frame.dat CRBNasLig.per-frame.dat)
receptor_prmtop=(2proteins.prmtop CRBN+protac.prmtop BRD4+protac.prmtop)
ligand_prmtop=(protac.prmtop BRD4.prmtop CRBN.prmtop)
trajectory="0-2.align.nc"
for ((i=0; i < 3 ; i++ )); do

mpirun -np 1 MMPBSA.py.MPI -O -i mmgbsa.in \
           -o  ${out[i]} \
           -eo ${out_frames[i]} \
           -cp ${complex_prmtop} \
           -rp ${receptor_prmtop[i]} \
           -lp ${ligand_prmtop[i]} \
            -y ${trajectory}

printf "\n${out[i]}\n" >> energy.dat
printf "\n\n\n\n\n${out_frames[i]}\n" >> energy.perframe.dat
cat ${out[i]} >> energy.dat
cat ${out_frames[i]} >> energy.perframe.dat

done
rm *en reference.frc *per-frame.dat mmgbsa.in
find -type l -delete


cd ../../

