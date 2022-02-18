#!/usr/bin/sh
#this script is used when the ligand is subject to protonation at the alkyl amine prior to the averaging.
for j in {1..18}; do
cd cluster$j ;


#pdbs=(combined*.pdb)

printf "parm 1/protac.pdb [1]\n" >cpptraj.average.in

for ((k=1; k < 40; k++)); do
if [ -d "$k" ]; then
printf "trajin $k/protac.pdb parm [1] #$k \n" >>cpptraj.average.in
fi
done
#printf "reference 1/2proteins.pdb parm [1]\nrms rms1 :1-149@CA reference \n" >>cpptraj.average.in
#printf "rms avgprotein.pdb\n" >>cpptraj.average.in

echo "average avgProtac.pdb pdb" >>cpptraj.average.in

cpptraj -i cpptraj.average.in

printf "parm 1/2proteins.pdb [1]\n" >cpptraj.average.in

for ((k=1; k < 40; k++)); do
if [ -d "$k" ]; then
printf "trajin $k/2proteins.pdb parm [1] #$k \n" >>cpptraj.average.in
fi
done


echo "average avgProtein.pdb pdb" >>cpptraj.average.in

cpptraj -i cpptraj.average.in

cat > convertmol2.cpptraj <<EOF
parm ../protac.prmtop

trajin avgProtac.pdb
trajout protac_gaff2.mol2
EOF
cpptraj convertmol2.cpptraj

mkdir ../rescore/cluster$j
mv avgProtac.pdb protac_gaff2.mol2 avgProtein.pdb ../rescore/cluster$j

cd ..
done
