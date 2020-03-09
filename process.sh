cat source_data/motifs/jaspar_genes2mat.txt | tr -d '\r' | sed -re 's/ +/\t/g' | sponge source_data/motifs/jaspar_genes2mat.txt
cat source_data/motifs/hocomoco_genes2mat.txt | tr -d '\r' | sed -re 's/ +/\t/g' | sponge source_data/motifs/hocomoco_genes2mat.txt

cat source_data/chipseq/remap_genes2exp.txt | tr -d '\r' | sed -re 's/ +/\t/g' | sponge source_data/chipseq/remap_genes2exp.txt
cat source_data/chipseq/remap_genes2exp.txt | tr -d '\r' | sed -re 's/ +/\t/g' | sponge source_data/chipseq/remap_genes2exp.txt

cat source_data/chipseq/hocomoco_remap_all_vs_all.txt | sed -re 's/^ +//' | sed -re 's/ +/\t/g' | sponge source_data/chipseq/hocomoco_remap_all_vs_all.txt
cat source_data/chipseq/jaspar_remap_all_vs_all.txt | sed -re 's/^ +//' | sed -re 's/ +/\t/g' | sponge source_data/chipseq/jaspar_remap_all_vs_all.txt
cat source_data/chipseq/hocomoco_remap_all_vs_all.txt <( tail -n+2 source_data/chipseq/jaspar_remap_all_vs_all.txt ) > source_data/chipseq/motifs_vs_remap.tsv

# SELEX aucs matrix was transposed compared to chipseq aucs matrix. We convert it to the same format
# Not idempotent!!! Do it once!!!
ruby transpose_selex_aucs_matrix.rb source_data/selex/hocomoco_jolma13_all_vs_all_roc10.txt | sponge source_data/selex/hocomoco_jolma13_all_vs_all_roc10.txt
ruby transpose_selex_aucs_matrix.rb source_data/selex/hocomoco_jolma13_all_vs_all_roc50.txt | sponge source_data/selex/hocomoco_jolma13_all_vs_all_roc50.txt
ruby transpose_selex_aucs_matrix.rb source_data/selex/jaspar_jolma13_all_vs_all_roc10.txt | sponge source_data/selex/jaspar_jolma13_all_vs_all_roc10.txt
ruby transpose_selex_aucs_matrix.rb source_data/selex/jaspar_jolma13_all_vs_all_roc50.txt | sponge source_data/selex/jaspar_jolma13_all_vs_all_roc50.txt

ruby transpose_selex_aucs_matrix.rb source_data/selex/hocomoco_jolma13ext_all_vs_all_roc10.txt | sponge source_data/selex/hocomoco_jolma13ext_all_vs_all_roc10.txt
ruby transpose_selex_aucs_matrix.rb source_data/selex/hocomoco_jolma13ext_all_vs_all_roc50.txt | sponge source_data/selex/hocomoco_jolma13ext_all_vs_all_roc50.txt
ruby transpose_selex_aucs_matrix.rb source_data/selex/jaspar_jolma13ext_roc10.txt | sponge source_data/selex/jaspar_jolma13ext_roc10.txt
ruby transpose_selex_aucs_matrix.rb source_data/selex/jaspar_jolma13ext_roc50.txt | sponge source_data/selex/jaspar_jolma13ext_roc50.txt

cat source_data/selex/hocomoco_jolma13_all_vs_all_roc10.txt <( tail -n+2 source_data/selex/jaspar_jolma13_all_vs_all_roc10.txt ) > source_data/selex/motifs_vs_selex10.tsv
cat source_data/selex/hocomoco_jolma13_all_vs_all_roc50.txt <( tail -n+2 source_data/selex/jaspar_jolma13_all_vs_all_roc50.txt ) > source_data/selex/motifs_vs_selex50.tsv

cat source_data/selex/hocomoco_jolma13ext_all_vs_all_roc10.txt <( tail -n+2 source_data/selex/jaspar_jolma13ext_roc10.txt ) > source_data/selex/motifs_vs_selex10.tsv
cat source_data/selex/hocomoco_jolma13ext_all_vs_all_roc50.txt <( tail -n+2 source_data/selex/jaspar_jolma13ext_roc50.txt ) > source_data/selex/motifs_vs_selex50.tsv

ruby download_jaspar_infos.rb
cat source_data/motifs/jaspar_infos/*.json > source_data/motifs/jaspar_infos.json
ruby download_jaspar_matrices.rb

mkdir source_data/motifs/all_pcms
cp source_data/motifs/hocomoco_pcm/* source_data/motifs/all_pcms/
cp source_data/motifs/jaspar_pcm/* source_data/motifs/all_pcms/
java -cp ape-3.0.2.jar ru.autosome.macroape.CollectDistanceMatrix source_data/motifs/all_pcms/ --pcm > distance_matrix.tsv
cat distance_matrix.tsv | ruby symmetrize_matrix.rb | sponge distance_matrix.tsv

ruby clusterize/clusterize.rb distance_matrix.tsv 0.95 --cluster-list clusters_dist_0.95.txt --mode distance
cut -f1 clusters_dist_0.95.txt | sort > source_data/motifs/representatives.txt

ruby collect_tf_annotation.rb

ruby annotate_cisbp.rb
ruby collect_motif_annotation.rb > source_data/annotation/motif_annotation_final.tsv
ruby extract_motif_features.rb source_data/annotation/motif_annotation_final.tsv > source_data/annotation/motif_annotation_final_features.tsv
( cat source_data/annotation/motif_annotation_final_features.tsv | head -1 | cuttab -f1,2,5,6,13,14,15;
  cat source_data/annotation/motif_annotation_final_features.tsv | tail -n+2 | cuttab -f1,2,5,6,13,14,15 | sort -u;
) > source_data/annotation/motif_features_uniq.tsv
ruby cleanup_annotation.rb

ruby motif_families.rb tf_class > motif_classes.tsv
ruby motif_families.rb tf_family > motif_families.tsv
ruby motif_families.rb cisbp_families > motif_cisbpfams.tsv

# best_motifs.rb was used to generate chipseq_best_motifs_with_classes.tsv and chipseq_best_motifs_with_families.tsv
ruby family_aggregator.rb chipseq_best_motifs_with_classes.tsv --drop-unknown-experiment --drop-unknown-motif > chipseq_best_motifs_with_classes_aggregated.tsv
# ruby family_aggregator.rb chipseq_best_motifs_with_families.tsv --drop-unknown-experiment --drop-unknown-motif > chipseq_best_motifs_with_families_aggregated.tsv

# See
# http://localhost:8080/alluvial-plot.html?tsv=chipseq_best_motifs_with_classes_aggregated.tsv

# C2H2 zinc finger factors{2.3.x}
cat chipseq_best_motifs_with_families.tsv | awk -F $'\t' -e '((NR==1) || ($6 ~ /\{2\.3\..+\}/) && ($2 != "unknown") && ($2 != "ambiguous")){print $0}'  > chipseq_best_motifs_families_C2H2-zinc-fingers.tsv

# Fork head / winged helix factors{3.3.x}
cat chipseq_best_motifs_with_families.tsv | awk -F $'\t' -e '((NR==1) || ($6 ~ /\{3\.3\..+\}/) && ($2 != "unknown") && ($2 != "ambiguous")){print $0}'  > chipseq_best_motifs_families_Forkhead.tsv

# Tryptophan cluster factors{3.5.x}
cat chipseq_best_motifs_with_families.tsv | awk -F $'\t' -e '((NR==1) || ($6 ~ /\{3\.5\..+\}/) && ($2 != "unknown") && ($2 != "ambiguous")){print $0}'  > chipseq_best_motifs_families_Tryptophan.tsv

./matrix_preparations.sh
./make_heatmaps.sh

################################
# Make heatmaps (old)
# mkdir -p results/chipseq/heatmap results/selex10/heatmap
# ruby mean_aucs_matrix.rb source_data/chipseq/motifs_vs_remap.tsv tf_family > results/chipseq/heatmap/heatmap_chipseq.tsv
# ruby mean_aucs_matrix.rb source_data/selex/motifs_vs_selex10.tsv tf_family > results/selex10/heatmap/heatmap_selex10.tsv

# ruby refine_heatmap.rb results/chipseq/heatmap/heatmap_chipseq.tsv tf_family > results/chipseq/heatmap/heatmap_chipseq.filtered.tsv
# ruby refine_heatmap.rb results/selex10/heatmap/heatmap_selex10.tsv tf_family > results/selex10/heatmap/heatmap_selex10.filtered.tsv

# ruby refine_heatmap.rb results/chipseq/heatmap/heatmap_chipseq.tsv tf_family --class-names > results/chipseq/heatmap/heatmap_chipseq.filtered.classnames.tsv
# ruby refine_heatmap.rb results/selex10/heatmap/heatmap_selex10.tsv tf_family --class-names > results/selex10/heatmap/heatmap_selex10.filtered.classnames.tsv

# python3.6 plot_heatmap.py results/chipseq/heatmap/heatmap_chipseq.filtered.tsv results/chipseq/heatmap/heatmap_chipseq.filtered.svg
# python3.6 plot_heatmap.py results/selex10/heatmap/heatmap_selex10.filtered.tsv results/selex10/heatmap/heatmap_selex10.filtered.svg

################################
# Find motifs which perform well

for AGGREGATION_LEVEL in tf_class tf_family cisbp_families; do
  ruby best_motifs.rb  source_data/final/remap_all_roc.txt $AGGREGATION_LEVEL --only-common-experiment-tfs > results/remap/best_motifs.${AGGREGATION_LEVEL}.tsv
  ruby best_motifs.rb  source_data/final/jy10_all_roc.txt  $AGGREGATION_LEVEL --only-common-experiment-tfs > results/jy10/best_motifs.${AGGREGATION_LEVEL}.tsv
  ruby best_motifs.rb  source_data/final/jy50_all_roc.txt  $AGGREGATION_LEVEL --only-common-experiment-tfs > results/jy50/best_motifs.${AGGREGATION_LEVEL}.tsv
  # ruby best_motifs.rb  source_data/chipseq/motifs_vs_remap.tsv  $AGGREGATION_LEVEL --only-common-experiment-tfs > results/chipseq/best_motifs.${AGGREGATION_LEVEL}.tsv
  # ruby best_motifs.rb  source_data/selex/motifs_vs_selex10.tsv  $AGGREGATION_LEVEL --only-common-experiment-tfs > results/selex10/best_motifs.${AGGREGATION_LEVEL}.tsv
  for AUC_THRESHOLD in 0.7  0.75  0.8; do
    ruby good_motifs.rb  source_data/final/remap_all_roc.txt  $AGGREGATION_LEVEL  $AUC_THRESHOLD --only-common-experiment-tfs > results/remap/good_motifs_${AUC_THRESHOLD}.${AGGREGATION_LEVEL}.tsv
    ruby good_motifs.rb  source_data/final/jy10_all_roc.txt  $AGGREGATION_LEVEL  $AUC_THRESHOLD --only-common-experiment-tfs > results/jy10/good_motifs_${AUC_THRESHOLD}.${AGGREGATION_LEVEL}.tsv
    ruby good_motifs.rb  source_data/final/jy50_all_roc.txt  $AGGREGATION_LEVEL  $AUC_THRESHOLD --only-common-experiment-tfs > results/jy50/good_motifs_${AUC_THRESHOLD}.${AGGREGATION_LEVEL}.tsv
    # ruby good_motifs.rb  source_data/chipseq/motifs_vs_remap.tsv  $AGGREGATION_LEVEL  $AUC_THRESHOLD --only-common-experiment-tfs > results/chipseq/good_motifs_${AUC_THRESHOLD}.${AGGREGATION_LEVEL}.tsv
    # ruby good_motifs.rb  source_data/selex/motifs_vs_selex10.tsv  $AGGREGATION_LEVEL  $AUC_THRESHOLD --only-common-experiment-tfs > results/selex10/good_motifs_${AUC_THRESHOLD}.${AGGREGATION_LEVEL}.tsv
  done
done

for AGGREGATION_LEVEL in tf_class tf_family cisbp_families; do
  ruby best_motifs.rb  source_data/final/uniprobe_all_cor.txt $AGGREGATION_LEVEL > results/uniprobe/best_motifs.${AGGREGATION_LEVEL}.tsv
  for COR_THRESHOLD in 0.2 0.3 0.4 0.5 0.6 0.7 0.75 0.8 0.85 0.9; do
    ruby good_motifs.rb  source_data/final/uniprobe_all_cor.txt  $AGGREGATION_LEVEL  $COR_THRESHOLD > results/uniprobe/good_motifs_${COR_THRESHOLD}.${AGGREGATION_LEVEL}.tsv
  done
done


# Open:
# python3 -m http.server 8080
#
# And see:
# http://localhost:8080/alluvial-plot.html?tsv=results/chipseq/best_motifs.classes.tsv
# http://localhost:8080/alluvial-plot.html?tsv=chipseq_best_motifs_families_Forkhead.tsv
# http://localhost:8080/alluvial-plot.html?tsv=chipseq_best_motifs_families_Tryptophan.tsv
# and
# http://localhost:8080/alluvial-plot-tfs.html?tsv=chipseq_best_motifs_families_C2H2-zinc-fingers.tsv
# http://localhost:8080/alluvial-plot-tfs.html?tsv=chipseq_best_motifs_families_Forkhead.tsv
# http://localhost:8080/alluvial-plot-tfs.html?tsv=chipseq_best_motifs_families_Tryptophan.tsv

for DATASET_TYPE in remap jy10 jy50; do
  ruby good_motifs.rb  source_data/final/${DATASET_TYPE}_all_roc.txt  cisbp_families  0.0 --only-common-experiment-tfs \
    | awktab -e '{print $4 "\t" $5 "\t" $1 "\t" $6 }' \
    | awktab -e '(NR==1 || $3 == $4){print $0}' \
    | ruby feature_quality_correlation.rb \
    > results/${DATASET_TYPE}/motif_tf_aucs.${DATASET_TYPE}.tsv

  ruby good_motifs.rb  source_data/final/${DATASET_TYPE}_all_roc.txt  cisbp_families  0.0 --only-common-experiment-tfs \
    | awktab -e '{print $4 "\t" $5 "\t" $1 "\t" $6 }' \
    | ruby feature_quality_correlation.rb \
    > results/${DATASET_TYPE}/motif_tf_aucs.all.${DATASET_TYPE}.tsv
done

# script needs manual fixes not to fail on missing values
ruby good_motifs.rb  source_data/final/uniprobe_all_cor.txt  cisbp_families  -100500.0 \
  | awktab -e '{print $4 "\t" $5 "\t" $1 "\t" $6 }' \
  | awktab -e '(NR==1 || $3 == $4){print $0}' \
  | ruby feature_quality_correlation.rb \
  > results/uniprobe/motif_tf_aucs.uniprobe.tsv

ruby good_motifs.rb  source_data/final/uniprobe_all_cor.txt  cisbp_families  -100500.0 \
  | awktab -e '{print $4 "\t" $5 "\t" $1 "\t" $6 }' \
  | ruby feature_quality_correlation.rb \
  > results/uniprobe/motif_tf_aucs.all.uniprobe.tsv


################################
ruby auc_infos_by_family.rb
ruby chipseq_selex_comparison.rb > chipseq_selex_comparison.tsv 
