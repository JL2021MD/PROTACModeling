##apply this script to the energy.per-frame.dat file to calculate the standard deviation.
##$1 is the file
filename=${1}
no_of_frames=100
cat -n $filename | grep -e Frame >& line_number

ln1=$(sed -n '1p' line_number | awk '{print $1+1}')
ln2=$(($ln1+$no_of_frames-1))

sed -n ""$ln1","$ln2"p" $filename | awk -F ',' '{print $NF }' > complex_total

ln1=$(sed -n '3p' line_number | awk '{print $1+1}')
ln2=$(($ln1+$no_of_frames-1))

sed -n ""$ln1","$ln2"p" $filename | awk -F ',' '{print $NF}' > protac_total

ln1=$(sed -n '6p' line_number | awk '{print $1+1}')
ln2=$(($ln1+$no_of_frames-1))

sed -n ""$ln1","$ln2"p" $filename | awk -F ',' '{print $NF}' > rec1

ln1=$(sed -n '7p' line_number | awk '{print $1+1}')
ln2=$(($ln1+$no_of_frames-1))

sed -n ""$ln1","$ln2"p" $filename | awk -F ',' '{print $NF}' > protein2_total

ln1=$(sed -n '10p' line_number | awk '{print $1+1}')
ln2=$(($ln1+$no_of_frames-1))

sed -n ""$ln1","$ln2"p" $filename | awk -F ',' '{print $NF}' > rec2

ln1=$(sed -n '11p' line_number | awk '{print $1+1}')
ln2=$(($ln1+$no_of_frames-1))

sed -n ""$ln1","$ln2"p" $filename | awk -F ',' '{print $NF}' > protein1_total

paste complex_total protac_total protein2_total protein1_total | awk '{print $1-$2-$3-$4}' > dH ; awk -v frames=$no_of_frames '{sum+=$1} END {print sum/frames}' dH > average_dH
ave=$(cat average_dH)
awk -v ave=$ave -v frames=$no_of_frames '{sum+=($1-ave)*($1-ave)} END {print sqrt(sum/(frames-1))}' dH > std_dH
std=$(cat std_dH)
echo $filename, "dH", $ave, $std

paste complex_total protac_total rec1 rec2 | awk '{print $1+$2-$3-$4}' > ddH ; awk -v frames=$no_of_frames '{sum+=$1} END {print sum/frames}' ddH > average_ddH
ave=$(cat average_ddH)
awk -v ave=$ave -v frames=$no_of_frames '{sum+=($1-ave)*($1-ave)} END {print sqrt(sum/(frames-1))}' ddH > std_ddH
std=$(cat std_ddH)
echo $filename, "ddH", $ave, $std

paste complex_total protein1_total protein2_total rec1 rec2 | awk '{print 2*$1-$2-$3-$4-$5 }' > Hcom ; awk -v frames=$no_of_frames '{sum+=$1} END {print sum/frames}' Hcom > average_Hcom
ave=$(cat average_Hcom)
awk -v ave=$ave -v frames=$no_of_frames '{sum+=($1-ave)*($1-ave)} END {print sqrt(sum/(frames-1))}' Hcom > std_Hcom
std=$(cat std_Hcom)
echo $filename, "Hcom", $ave, $std
