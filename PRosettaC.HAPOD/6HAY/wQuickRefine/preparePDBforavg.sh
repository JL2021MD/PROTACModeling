#this script is used when the ligand needs to be reprotonated. Because the alkyl amine may not be properly protonated with a +1 net charge. The script extracts the seperate ligand from each PRosettaC pose, protonates it, and generates its average pose.
for i in {1..16}; do
cd cluster$i
pdbs=combined_*
cat >extractligand.py<<EOF
import os
import sys
import re
#print("$i")
with open('ligand.pdb', 'w') as f, open('2proteinsWithH.pdb', 'w') as p, open(sys.argv[1],"r") as file:
	lines = file.readlines()
	for index, line in enumerate(lines):
		if "PTC" in line:
			print(lines[index].strip(), file=f)
			
		elif "ATOM" and " A " in line:
			print(lines[index].strip(), file=p)
	for index, line in enumerate(lines):
		if "ATOM" and " B " in line:
			print(lines[index].strip(), file=p)
		
EOF
k=1
for j in ${pdbs[@]}; do
python3 extractligand.py $j
mkdir $k
mv ligand.pdb 2proteinsWithH.pdb $k/

cd $k
cat >delHprotein.com<<EOF

sel H
del sel
write #0 2proteins.pdb

EOF

cat >addHligand.com<<EOF

addh hbond false
write #0 protac.pdb

EOF


chimera --nogui 2proteinsWithH.pdb delHprotein.com
chimera --nogui ligand.pdb addHligand.com
k=$((k+1))
rm *.com *.py
cd ..
done
cd ..
done
exit

cat > convertmol2.cpptraj <<EOF
parm ../protac.prmtop

trajin protac.pdb
trajout protac_gaff2.mol2
EOF
cpptraj convertmol2.cpptraj
done
exit
mkdir ../rescore/cluster$i

mv protac.pdb protac_gaff2.mol2 2proteins.pdb ../rescore/cluster$i/
#cp rep.c0.pdb ../rescore/cluster$i/

cd ..


done
