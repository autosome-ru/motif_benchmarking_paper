cat source_data/motifs/jaspar_genes2mat.txt | tr -d '\r' | sed -re 's/ +/\t/g' | sponge source_data/motifs/jaspar_genes2mat.txt
cat source_data/motifs/hocomoco_genes2mat.txt | tr -d '\r' | sed -re 's/ +/\t/g' | sponge source_data/motifs/hocomoco_genes2mat.txt

cat source_data/chipseq/remap_genes2exp.txt | tr -d '\r' | sed -re 's/ +/\t/g' | sponge source_data/chipseq/remap_genes2exp.txt

cat source_data/chipseq/hocomoco_remap_all_vs_all.txt | sed -re 's/^ +//' | sed -re 's/ +/\t/g' | sponge source_data/chipseq/hocomoco_remap_all_vs_all.txt
cat source_data/chipseq/jaspar_remap_all_vs_all.txt | sed -re 's/^ +//' | sed -re 's/ +/\t/g' | sponge source_data/chipseq/jaspar_remap_all_vs_all.txt
cat source_data/chipseq/hocomoco_remap_all_vs_all.txt <( tail -n+2 source_data/chipseq/jaspar_remap_all_vs_all.txt ) > source_data/chipseq/motifs_vs_remap.tsv

ruby collect_tf_annotation.rb
ruby motif_families.rb 2 > motif_classes.tsv
ruby motif_families.rb 3 > motif_families.tsv
ruby best_motifs.rb  source_data/chipseq/motifs_vs_remap.tsv  2 > chipseq_best_motifs_with_classes.tsv
ruby best_motifs.rb  source_data/chipseq/motifs_vs_remap.tsv  3 > chipseq_best_motifs_with_families.tsv

ruby family_aggregator.rb chipseq_best_motifs_with_classes.tsv --drop-unknown-experiment --drop-unknown-motif > chipseq_best_motifs_with_classes_aggregated.tsv
# ruby family_aggregator.rb chipseq_best_motifs_with_families.tsv --drop-unknown-experiment --drop-unknown-motif > chipseq_best_motifs_with_families_aggregated.tsv

# C2H2 zinc finger factors{2.3.x}
cat chipseq_best_motifs_with_families.tsv | awk -F $'\t' -e '((NR==1) || ($6 ~ /\{2\.3\..+\}/) && ($2 != "unknown") && ($2 != "ambiguous")){print $0}'  > chipseq_best_motifs_families_C2H2-zinc-fingers.tsv

# Fork head / winged helix factors{3.3.x}
cat chipseq_best_motifs_with_families.tsv | awk -F $'\t' -e '((NR==1) || ($6 ~ /\{3\.3\..+\}/) && ($2 != "unknown") && ($2 != "ambiguous")){print $0}'  > chipseq_best_motifs_families_Forkhead.tsv

# Tryptophan cluster factors{3.5.x}
cat chipseq_best_motifs_with_families.tsv | awk -F $'\t' -e '((NR==1) || ($6 ~ /\{3\.5\..+\}/) && ($2 != "unknown") && ($2 != "ambiguous")){print $0}'  > chipseq_best_motifs_families_Tryptophan.tsv
