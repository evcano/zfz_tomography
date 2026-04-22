MCTOMODIR="/home/valeroe/packages/MCTomo/utils"
INVDIR="/data/valeroe/red_sea_obs_2/3d_inversion/final_chains_2million"

# iterate over chains
for i in {01..20}
do
    CDIR="chain_"$i
    echo $CDIR
    cd $INVDIR/$CDIR/Results
    cp $MCTOMODIR/samp2txt .
    pwd
    echo "running samp2txt "
    ./samp2txt
    echo "samp2txt finished"
    echo "returning to 3d inv dir"
    cd $INVDIR
done
