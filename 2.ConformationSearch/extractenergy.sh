##run this script one folder above the trajectories folder. Rename the folders that are not by a simple number only
printf "trajectory\ncomplex\n2proteins\nprotac\ncomplex\nrec1\nprotein1\ncomplex\nrec2\nprotein2" > energyPMEMD0 
for i in PMEMD* ; do

printf "$i \n" >> energy$i

awk '/TOTAL                  / { print $2 }' $i/energy.results.dat >> energy$i

#printf "\n" >> energy_extracted.csv

done

paste -d "," energyPMEMD{1..30} > energycombined.csv
rm energyPMEMD*
