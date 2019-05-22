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

cat source_data/selex/hocomoco_jolma13_all_vs_all_roc10.txt <( tail -n+2 source_data/selex/jaspar_jolma13_all_vs_all_roc10.txt ) > source_data/selex/motifs_vs_selex10.tsv
cat source_data/selex/hocomoco_jolma13_all_vs_all_roc50.txt <( tail -n+2 source_data/selex/jaspar_jolma13_all_vs_all_roc50.txt ) > source_data/selex/motifs_vs_selex50.tsv

ruby download_jaspar_matrices.rb

mkdir source_data/motifs/all_pcms
cp source_data/motifs/hocomoco_pcm/* source_data/motifs/all_pcms/
cp source_data/motifs/jaspar_pcm/* source_data/motifs/all_pcms/
java -cp ape-3.0.2.jar ru.autosome.macroape.CollectDistanceMatrix source_data/motifs/all_pcms/ --pcm > distance_matrix.tsv
cat distance_matrix.tsv | ruby symmetrize_matrix.rb | sponge distance_matrix.tsv

ruby clusterize/clusterize.rb distance_matrix.tsv 0.95 --cluster-list clusters_dist_0.95.txt --mode distance
cut -f1 clusters_dist_0.95.txt | sort > source_data/motifs/representatives.txt

ruby collect_tf_annotation.rb
ruby motif_families.rb 2 > motif_classes.tsv
ruby motif_families.rb 3 > motif_families.tsv
ruby best_motifs.rb  source_data/chipseq/motifs_vs_remap.tsv  2 > chipseq_best_motifs_with_classes.tsv
ruby best_motifs.rb  source_data/chipseq/motifs_vs_remap.tsv  3 > chipseq_best_motifs_with_families.tsv

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

# See
# http://localhost:8080/alluvial-plot.html?tsv=chipseq_best_motifs_families_C2H2-zinc-fingers.tsv
# http://localhost:8080/alluvial-plot.html?tsv=chipseq_best_motifs_families_Forkhead.tsv
# http://localhost:8080/alluvial-plot.html?tsv=chipseq_best_motifs_families_Tryptophan.tsv
# and
# http://localhost:8080/alluvial-plot-tfs.html?tsv=chipseq_best_motifs_families_C2H2-zinc-fingers.tsv
# http://localhost:8080/alluvial-plot-tfs.html?tsv=chipseq_best_motifs_families_Forkhead.tsv
# http://localhost:8080/alluvial-plot-tfs.html?tsv=chipseq_best_motifs_families_Tryptophan.tsv

ruby auc_infos_by_family.rb
