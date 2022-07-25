for i in {1..5}; do
cd $i
rm add_charge_frcmod.sh
ln -s ../tleap_HMR.sh
bash tleap_HMR.sh
mkdir Refine{1..5}
cd ..
done
