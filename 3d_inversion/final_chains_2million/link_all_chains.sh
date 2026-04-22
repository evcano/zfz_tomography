INVDIR="/data/valeroe/red_sea_obs/3d_inversion/final_chains_2million"
OUTDIR="/data/valeroe/red_sea_obs/3d_inversion/final_chains_2million/all_chains_TEST"

for i in {01..20}
do

j=$(expr $i + 0)

# ini
ln -s $INVDIR"/chain_"$i"/Results/ini/InitialSample_1.dat" $OUTDIR"/Results/ini/InitialSample_"$j".dat"
ln -s $INVDIR"/chain_"$i"/Results/ini/InitialSample_1.txt" $OUTDIR"/Results/ini/InitialSample_"$j".txt"
ln -s $INVDIR"/chain_"$i"/Results/ini/InitialSigma_1.dat"  $OUTDIR"/Results/ini/InitialSigma_"$j".dat"

# last run
ln -s $INVDIR"/chain_"$i"/Results/last_run/last_average_1.dat"  $OUTDIR"/Results/last_run/last_average_"$j".dat"
ln -s $INVDIR"/chain_"$i"/Results/last_run/last_info_1.dat"     $OUTDIR"/Results/last_run/last_info_"$j".dat"
ln -s $INVDIR"/chain_"$i"/Results/last_run/last_sources_1.dat"  $OUTDIR"/Results/last_run/last_sources_"$j".dat"
ln -s $INVDIR"/chain_"$i"/Results/last_run/last_var_1.dat"      $OUTDIR"/Results/last_run/last_var_"$j".dat"
ln -s $INVDIR"/chain_"$i"/Results/last_run/last_vertices_1.dat" $OUTDIR"/Results/last_run/last_vertices_"$j".dat"
ln -s $INVDIR"/chain_"$i"/Results/last_run/last_vertices_1.txt" $OUTDIR"/Results/last_run/last_vertices_"$j".txt"

#
ln -s $INVDIR"/chain_"$i"/Results/likelihood_1.dat" $OUTDIR"/Results/likelihood_"$j".dat"
ln -s $INVDIR"/chain_"$i"/Results/ncells_1.dat"     $OUTDIR"/Results/ncells_"$j".dat"
ln -s $INVDIR"/chain_"$i"/Results/samples_1.out"    $OUTDIR"/Results/samples_"$j".out"

done
