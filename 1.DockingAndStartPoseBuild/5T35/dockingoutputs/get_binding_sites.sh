#!/bin/bash
# this script exctracts the raw coordinates of selected atoms from each PDB file. 

printf "name,lig_atom1_x,y,z,lig_atom2_x,y,z,lig_atom3_x,y,z,rec_atom_1_x,y,z,rec_atom_2_x,y,z,rec_atom3_x,y,z \n" >>distance.csv
for inpdb in *.pdb; do


# match these strings with grep to extract atom lines for binding site
#VHL
lig1="CA  ASN C  67"
lig2="CA  PHE C  91"
lig3="CA  ARG C  69"

#BRD4
rec1="CA  ASN A 433"
rec2="CA  LEU A 385"
rec3="CA  TYR A 432"

# grep will find lines that match strings, awk will extract the x, y, z columns from each line
ligX=`cat $inpdb | grep -e "${lig1}" -e "${lig2}" -e "${lig3}" | awk '{print $7","$8","$9","}'`
recX=`cat $inpdb | grep -e "${rec1}" -e "${rec2}" -e "${rec3}" | awk '{print $7","$8","$9","}'`

# print to file
 
echo ${inpdb}","${ligX} ${recX} >> distance.csv
done

#The distance.csv file was opened in a spreadsheets editor using comma as the deliminator, and the average coordinates was calculated within the spreadsheets, followed by calculating of the distance.
