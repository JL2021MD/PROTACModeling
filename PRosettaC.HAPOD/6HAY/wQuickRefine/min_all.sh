for i in {1..16} ; do
cd cluster$i
mkdir min
cd min
bash ../../min.sh
cd ../..
done
