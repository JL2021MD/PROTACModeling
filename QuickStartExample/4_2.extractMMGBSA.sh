## MMGBSA is used as a prescore in the original publication, mainly for poses(distinct states) found from long MD relaxations. However, if poses directly from ternary docking is used as input, as in the current example, the performance of MMGBSA, even for prescoring may drop significantly. Based on single-point enthalpy calculations after minimization of the docked pose, the near-native pose 11_11 does not score favorable at all. It is therefore suggested to not use MMGBSA at all to score ternary docking poses. By comparison, long MD relaxations are expected to remove false positives and have the binding configuration converge onto distinct states of better local minima. Those states (and not a single pose) are better suited for MMGBSA prescoring. 

## The following script extracts the calculated single-point MMGBSA energy calculated for the 5 starting poses after minimization.

printf "Pose\ncomplex\n2proteins\nprotac\ncomplex\nrec1\nprotein1\ncomplex\nrec2\nprotein2" > energy0 
for i in {1..5} ; do

printf "$i \n" >> energy$i

awk '/TOTAL                  / { print $2 }' $i/formin/energy.dat >> energy$i

done

if paste -d "," energy{0..5} > energycombined.csv  ; then
printf "Energies combined into >> energycombined.csv \n"
fi
rm energy{0..5}
