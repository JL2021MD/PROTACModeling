#!/bin/bash
## Adding AM1-BCC charges to the PROTAC molecules may take 1-2 hours, because of the molecule being pretty large. For the example run, to save time one can use the already generated output files here. Set demo=0 to run yourself.

demo=0
if [ $demo -eq 0 ]
then
for i in {1..4}; do
cd $i
cat >add_charge_frcmod.sh <<EOF

antechamber -i protac.pdb -fi pdb -o protac_gaff2.mol2 -fo mol2 -c bcc -nc 0 -at gaff2 

rm *.AC *.INF *.AC0 sqm.pdb sqm.in sqm.out
parmchk2 -i protac_gaff2.mol2 -f mol2 -s 2 -o ligand.frcmod
printf "Finished preparing Pose $i ligand.\n\n"
EOF
bash add_charge_frcmod.sh &
cd ..
done

##They will run in parallel in the background to save time. The last loop is written seperately so the the script can end.
cd 5
cat >add_charge_frcmod.sh <<EOF

antechamber -i protac.pdb -fi pdb -o protac_gaff2.mol2 -fo mol2 -c bcc -nc 0 -at gaff2 

rm *.AC *.INF *.AC0 sqm.pdb sqm.in sqm.out
parmchk2 -i protac_gaff2.mol2 -f mol2 -s 2 -o ligand.frcmod

printf "Finished preparing Pose 5 ligand.\n\n"
EOF
bash add_charge_frcmod.sh
cd ..

echo "Wait for all background processes (sqm) to finish and all (5) poses to print the finished messages."
# -nc 0 is the net charge of the ligand being 0 in our case.

else
for i in {1..5}; do
ln -s ../Other/ForQuickDemo/$i/{protac_gaff2.mol2,ligand.frcmod} $i/ -v
done

fi
exit

## to undo, run the following
for i in {1..5}; do
rm -v $i/{protac_gaff2.mol2,ligand.frcmod}
done
