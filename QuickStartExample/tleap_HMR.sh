current=`pwd`
mkdir $current/Tleap
cd  $current/Tleap

cat > tleap.build.in <<EOF
#### tleap -f tleap.build.in


source leaprc.protein.ff19SB
source leaprc.gaff2

source leaprc.water.opc
loadamberparams frcmod.ions1lm_126_iod_opc

set default PBradii mbondi2

rec=loadpdb ../2proteins.pdb

###Load Ligand frcmod/prep

loadamberparams ../ligand.frcmod
lig=loadmol2 ../protac_gaff2.mol2

gascomplex= combine {rec lig}

savepdb gascomplex complex.gas.pdb
saveamberparm gascomplex complex.gas.prmtop complex.gas.rst7

solvcomplex= combine {rec lig}
solvateoct solvcomplex OPCBOX 10.0

###Neutralize system (it will add either Na or Cl depending on net charge)
addions solvcomplex Cl- 0
addions solvcomplex Na+ 0

###Write solvated pdb file
savepdb solvcomplex complex.opc.pdb

charge solvcomplex

###Write Solvated topology and coord file
saveamberparm solvcomplex complex.opc.prmtop complex.opc.rst7

quit

EOF

cat > parmed1.in<<EOF
HMassRepartition
outparm complex.HMR.opc.prmtop
EOF
cat > parmed2.in<<EOF
HMassRepartition
outparm complex.HMR.gas.prmtop
EOF
tleap -f tleap.build.in
parmed -p complex.opc.prmtop -i parmed1.in

cat > stripPrmtop.temp <<EOF

parm complex.gas.prmtop [1]
parm complex.gas.prmtop [2]
parm complex.gas.prmtop [3]
parm complex.gas.prmtop [4]
parm complex.gas.prmtop [5]
parm complex.gas.prmtop [6]

parmstrip :1-127:512 parm [1]
parmwrite out CRBN.prmtop parm [1]

parmstrip :1-511 parm [2]
parmwrite out protac.prmtop parm [2]

parmstrip :128-512 parm [3]
parmwrite out BRD4.prmtop parm [3]

parmstrip :512 parm [4]
parmwrite out 2proteins.prmtop parm [4]

parmstrip :1-127 parm [5]
parmwrite out CRBN+protac.prmtop parm [5]

parmstrip :128-511 parm [6]
parmwrite out BRD4+protac.prmtop parm [6]

EOF
cat > zinc.in <<EOF
parm complex.gas.prmtop
reference complex.gas.rst7
rst :475@SG :407@SG reference offset 0.0 rk2 0.0 rk3 5.0 out zinc.rst
rst :478@SG :410@SG reference offset 0.0 rk2 0.0 rk3 5.0 out zinc.rst

EOF

cpptraj -i stripPrmtop.temp -o cpptraj.log.temp
cpptraj -i zinc.in -o cpptraj.log.temp

