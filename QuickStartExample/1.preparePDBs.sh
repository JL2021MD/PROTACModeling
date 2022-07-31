##If there are protonation states for amino acid side chains different than the common ones, they should be edited first (e.g. HID). In this example run the His378 in CRBN has been edited into HID.

##The zinc atom was deleted during the ternary docking process. It can be added back by alignment. However it has minimum effect on the HAPOD scoring, and since adding distance restraints between the zinc finger sulfurs will have the same effect, so for simplicity the zinc is not added back in this example.

pdbs=cluster*.pdb
##This splits the input PDB into the small molecule and 2 proteins.
cat >extractligand.py<<EOF
import os
import sys
import re

with open('ligand.pdb', 'w') as f, open('2proteinsWithH.pdb', 'w') as p, open(sys.argv[1],"r") as file:
	lines = file.readlines()
	for index, line in enumerate(lines):
		if "UNL" in line:
			print(lines[index].strip(), file=f)
		elif "CONECT" in line:
			print(lines[index].strip(), file=f)
## "UNL" is your protac residue ID, change this accordingly
		if "ATOM" and " B " in line:
			print(lines[index].strip(), file=p)
	for index, line in enumerate(lines):
## writing chain A and B seperately enables the possiblity of switching their order if needed
		if "ATOM" and " A " in line:
			print(lines[index].strip(), file=p)
		
EOF
k=1
for j in ${pdbs[@]}; do
python3 extractligand.py $j
mkdir $k
mv ligand.pdb 2proteinsWithH.pdb $k/

cd $k
ln -s ../$j .
antechamber -i ligand.pdb -fi pdb -o ligand2.pdb -fo pdb
#this gives atoms unique names, needed for selecting atoms in Chimera (next script) to correct the atom type.

#the selected atoms here pertain to the ones that need atom type correction in the demo, after being already inspected in Chimera. For your own molecules, run this 1.preparePDBs.sh once without the 2 lines of "setattr" below, and inspect to see which atoms need correction, and add to the "setattr" lines below. Delete the 1-xx folders, and then rerun 1.preparePDBs.sh. You may need to do this several times iteratively. Most errors will be caught in the 2.prepareForRefine.sh step.
cat >addH.com<<EOF
setattr a idatmType C3 :UNL@C,C2,C3,C4,C5,C6,C7,C8,C9
setattr a idatmType Npl :UNL@N1,N
sel H
del sel
addh
write #1 protac.pdb
sel H
del sel
write #0 2proteins.pdb
EOF

## UCSF Chimera is used here for adding or deleting hydrogens. It is also a decent and free program to visualize structures and easily generate beautiful images. 

chimera --nogui 2proteinsWithH.pdb ligand2.pdb addH.com
k=$((k+1))
rm *.com

#clean up
rm ligand2.pdb ligand.pdb
cd ..
done
rm *.py

