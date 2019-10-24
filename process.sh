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

ruby motif_families.rb 2 > motif_classes.tsv
ruby motif_families.rb 3 > motif_families.tsv

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
ruby best_motifs.rb  source_data/chipseq/motifs_vs_remap.tsv  2  > results/chipseq/best_motifs.classes.tsv
ruby best_motifs.rb  source_data/chipseq/motifs_vs_remap.tsv  3  > results/chipseq/best_motifs.families.tsv
ruby best_motifs.rb  source_data/selex/motifs_vs_selex10.tsv  2  > results/selex10/best_motifs.classes.tsv
ruby best_motifs.rb  source_data/selex/motifs_vs_selex10.tsv  3  > results/selex10/best_motifs.families.tsv

ruby good_motifs.rb  source_data/chipseq/motifs_vs_remap.tsv  2  0.7  > results/chipseq/good_motifs_0.7.classes.tsv
ruby good_motifs.rb  source_data/chipseq/motifs_vs_remap.tsv  2  0.75  > results/chipseq/good_motifs_0.75.classes.tsv
ruby good_motifs.rb  source_data/chipseq/motifs_vs_remap.tsv  2  0.8  > results/chipseq/good_motifs_0.8.classes.tsv

ruby good_motifs.rb  source_data/chipseq/motifs_vs_remap.tsv  3  0.7  > results/chipseq/good_motifs_0.7.families.tsv
ruby good_motifs.rb  source_data/chipseq/motifs_vs_remap.tsv  3  0.75  > results/chipseq/good_motifs_0.75.families.tsv
ruby good_motifs.rb  source_data/chipseq/motifs_vs_remap.tsv  3  0.8  > results/chipseq/good_motifs_0.8.families.tsv

ruby good_motifs.rb  source_data/selex/motifs_vs_selex10.tsv  2  0.7  > results/selex10/good_motifs_0.7.classes.tsv
ruby good_motifs.rb  source_data/selex/motifs_vs_selex10.tsv  2  0.75  > results/selex10/good_motifs_0.75.classes.tsv
ruby good_motifs.rb  source_data/selex/motifs_vs_selex10.tsv  2  0.8  > results/selex10/good_motifs_0.8.classes.tsv

ruby good_motifs.rb  source_data/selex/motifs_vs_selex10.tsv  3  0.7  > results/selex10/good_motifs_0.7.families.tsv
ruby good_motifs.rb  source_data/selex/motifs_vs_selex10.tsv  3  0.75  > results/selex10/good_motifs_0.75.families.tsv
ruby good_motifs.rb  source_data/selex/motifs_vs_selex10.tsv  3  0.8  > results/selex10/good_motifs_0.8.families.tsv

# See
# http://localhost:8080/alluvial-plot.html?tsv=results/chipseq/best_motifs.classes.tsv
# http://localhost:8080/alluvial-plot.html?tsv=chipseq_best_motifs_families_Forkhead.tsv
# http://localhost:8080/alluvial-plot.html?tsv=chipseq_best_motifs_families_Tryptophan.tsv
# and
# http://localhost:8080/alluvial-plot-tfs.html?tsv=chipseq_best_motifs_families_C2H2-zinc-fingers.tsv
# http://localhost:8080/alluvial-plot-tfs.html?tsv=chipseq_best_motifs_families_Forkhead.tsv
# http://localhost:8080/alluvial-plot-tfs.html?tsv=chipseq_best_motifs_families_Tryptophan.tsv

################################
ruby auc_infos_by_family.rb
ruby chipseq_selex_comparison.rb > chipseq_selex_comparison.tsv 
