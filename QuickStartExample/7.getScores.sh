set -e

for i in {1..5} ; do
cd $i/HAPOD

for output in Occupancy.out Temperature.out; do
if [ -e $output ]; then
rm $output
fi
done
for i in PMEMD* ; do

printf "$i, " >> Occupancy.out
printf "$i, " >> Temperature.out

awk '/Occupancy/ { print $3 }'  $i/HAPOD.out >> Occupancy.out
awk '/Temperature/ { print $3 }' $i/HAPOD.out >> Temperature.out

done

current=`pwd`

for output in Occupancy.out Temperature.out; do
printf "\nData from directory $current. \n"
cat > stats.py <<EOF
import sys
input=open(sys.argv[1], "r")
lines=input.readlines()
scores=[]
for index,line in enumerate(lines):
    value=float(line.split()[1])
    scores.append(value)

l=len(scores)
average=sum(scores)/l
print("Average",sys.argv[1].rsplit(".")[0],"score of", l, "runs is:", round(average,1))
var=sum([(value-average)**2 for value in scores]) /l

dev=var**0.5
print("Standard deviation: ", round(dev,1))

err=dev/l**0.5
print("Standard error: ", round(err,1))

EOF
python3 stats.py $output
done
cd ../../
done

printf "\nThank you for trying HAPOD scoring, if it is helpful in any way, please cite: \nLiao J, Nie X, Unarta IC, Ericksen SS, Tang W. In Silico Modeling and Scoring of PROTAC-Mediated Ternary Complex Poses. J Med Chem. 2022 Apr 28;65(8):6116-6132. doi: 10.1021/acs.jmedchem.1c02155.\n\n"

