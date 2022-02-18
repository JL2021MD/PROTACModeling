#This helps combine all the MMGBSA energies of each averaged and minimized frame of from each PRosettaC cluster.
for i in {1..20} ; do

printf "$i \n" >> energy_extracted.csv

awk '/TOTAL                  / { print $2 }' cluster$i/min/energy.min.dat >> energy_extracted.csv

printf "\n" >> energy_extracted.csv

done
