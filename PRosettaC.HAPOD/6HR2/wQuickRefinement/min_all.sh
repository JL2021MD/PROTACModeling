for i in {1..18} ; do
cd cluster$i/min
bash min.sh
cd ../..
done
