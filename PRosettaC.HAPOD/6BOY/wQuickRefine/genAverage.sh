#This is run in each folder of each PRosettaC cluster.
pdbs=(*.pdb)

k=1
for i in *pdb ; do
printf "parm $i [$k]\n" >>cpptraj.average.in
printf "trajin $i parm [$k] #$k\n" >>cpptraj.average.in
k=$((k+1))
done
printf "trajout average.pdb pdb" >>cpptraj.average.in
cpptraj cpptraj.average.in
