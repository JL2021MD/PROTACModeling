##This optional script helps users quickly get a sense of how different each input pose for HAPOD rescore is from the native pose (reference file).

##The reference file is made by alignments (in Chimera) to the crystal pose, from one of the output poses. The crystal pose shouldn't be directly used as the reference because it will have different lengths and missing chains causing inconsistency.

##Most poses using Integrative PROTAC-Model for ternary docking will be from RosettaDock, but a few will be from the Frodock stage. In our example here, 24_24 is from the Frodock stage. This results in 24_24 being a few atoms less than the rest. cpptraj can use a slightly shorter protein as the parm file, but not the other way around. So cluster_24_24.noH.pdb is the parm file. Hydrogens are removed the futher minimize the effects of atom naming mismatches. Alternatively, omit 24_24 from the RMSD measurements and choose any other of the 4 PDBs as parm.

##Please note: The terms RosettaDock, Frodock are tools used in the ternary docking. For running this quick start example you do not need to get involved in those steps, however, since its output is already included here.

cd .. 
ln -s Other/reference.pdb
for i in *pdb ; do
reduce -Trim $i > ${i%.pdb}.noH.pdb
done
pdbs=(*.noH.pdb)
printf "parm ${pdbs[3]} [1]\n" > cpptraj.RMSD.in
k=1
for i in *.noH.pdb ; do
printf "trajin $i parm [1] #$k\n" >>cpptraj.RMSD.in
k=$((k+1))
done

printf "reference reference.noH.pdb  \n" >>cpptraj.RMSD.in
printf "rmsd CRBN reference ":1-384@CA" out RMSDresults.out nofit\n" >>cpptraj.RMSD.in
printf "rmsd BRD4 reference ":385-511@CA" out RMSDresults.out nofit\n" >>cpptraj.RMSD.in
printf "rmsd bestfit(entire) reference ":1-511@CA" out RMSDresults.out nofit\n" >>cpptraj.RMSD.in

cpptraj cpptraj.RMSD.in
rm *.noH.pdb cpptraj.RMSD.in
mv RMSDresults.out ./Other
find -type l -delete
