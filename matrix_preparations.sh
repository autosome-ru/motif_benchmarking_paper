mkdir -p source_data/final/

#############

ruby transpose_selex_aucs_matrix.rb ./ccg.epfl.ch/pwmtools/benchmarking/remap_hocomoco_roc.txt > source_data/final/remap_hocomoco_roc.txt
ruby transpose_selex_aucs_matrix.rb ./ccg.epfl.ch/pwmtools/benchmarking/remap_jaspar_roc.txt   > source_data/final/remap_jaspar_roc.txt
ruby transpose_selex_aucs_matrix.rb ./ccg.epfl.ch/pwmtools/benchmarking/remap_cisbp_roc.txt    > source_data/final/remap_cisbp_roc.txt


cat source_data/final/remap_hocomoco_roc.txt | head -1    > source_data/final/remap_all_roc.txt
cat source_data/final/remap_hocomoco_roc.txt | tail -n+2 >> source_data/final/remap_all_roc.txt
cat source_data/final/remap_jaspar_roc.txt   | tail -n+2 >> source_data/final/remap_all_roc.txt
cat source_data/final/remap_cisbp_roc.txt    | tail -n+2 >> source_data/final/remap_all_roc.txt

#############
for JY_TOP in 10 50; do
    ruby transpose_selex_aucs_matrix.rb ./ccg.epfl.ch/pwmtools/benchmarking/jolma_yang_hocomoco_roc${JY_TOP}.txt > source_data/final/jy${JY_TOP}_hocomoco_roc.txt
    ruby transpose_selex_aucs_matrix.rb ./ccg.epfl.ch/pwmtools/benchmarking/jolma_yang_jaspar_roc${JY_TOP}.txt   > source_data/final/jy${JY_TOP}_jaspar_roc.txt
    ruby transpose_selex_aucs_matrix.rb ./ccg.epfl.ch/pwmtools/benchmarking/jolma_yang_cisbp_roc${JY_TOP}.txt    > source_data/final/jy${JY_TOP}_cisbp_roc.txt

    cat source_data/final/jy${JY_TOP}_hocomoco_roc.txt | head -1    > source_data/final/jy${JY_TOP}_all_roc.txt
    cat source_data/final/jy${JY_TOP}_hocomoco_roc.txt | tail -n+2 >> source_data/final/jy${JY_TOP}_all_roc.txt
    cat source_data/final/jy${JY_TOP}_jaspar_roc.txt   | tail -n+2 >> source_data/final/jy${JY_TOP}_all_roc.txt
    cat source_data/final/jy${JY_TOP}_cisbp_roc.txt    | tail -n+2 >> source_data/final/jy${JY_TOP}_all_roc.txt
done

#############

cat source_data/final/LOG_correlations_hs-mm_manually_selected.tsv | head -1 | sed -re 's/.txt\b//g'  > source_data/final/uniprobe_all_cor.txt
cat source_data/final/LOG_correlations_hs-mm_manually_selected.tsv | ruby sort_columns.rb | tail -n+2 >> source_data/final/uniprobe_all_cor.txt
cat source_data/final/LOG_correlations_cisbp_hs-mm_manually_selected.tsv | ruby sort_columns.rb | tail -n+2 >> source_data/final/uniprobe_all_cor.txt
