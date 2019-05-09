cat source_data/motifs/jaspar_genes2mat.txt | tr -d '\r' | sed -re 's/ +/\t/g' | sponge source_data/motifs/jaspar_genes2mat.txt
cat source_data/motifs/hocomoco_genes2mat.txt | tr -d '\r' | sed -re 's/ +/\t/g' | sponge source_data/motifs/hocomoco_genes2mat.txt

cat source_data/chipseq/remap_genes2exp.txt | tr -d '\r' | sed -re 's/ +/\t/g' | sponge source_data/chipseq/remap_genes2exp.txt

cat source_data/chipseq/hocomoco_remap_all_vs_all.txt | sed -re 's/^ +//' | sed -re 's/ +/\t/g' | sponge source_data/chipseq/hocomoco_remap_all_vs_all.txt
cat source_data/chipseq/jaspar_remap_all_vs_all.txt | sed -re 's/^ +//' | sed -re 's/ +/\t/g' | sponge source_data/chipseq/jaspar_remap_all_vs_all.txt
cat source_data/chipseq/hocomoco_remap_all_vs_all.txt <( tail -n+2 source_data/chipseq/jaspar_remap_all_vs_all.txt ) > source_data/chipseq/motifs_vs_remap.tsv

ruby collect_tf_annotation.rb
ruby motif_families.rb 2 > motif_classes.tsv
ruby motif_families.rb 3 > motif_families.tsv
ruby best_motifs.rb  source_data/chipseq/motifs_vs_remap.tsv  source_data/chipseq/remap_genes2exp.txt  2 > chipseq_best_motifs_with_classes.tsv
ruby best_motifs.rb  source_data/chipseq/motifs_vs_remap.tsv  source_data/chipseq/remap_genes2exp.txt  3 > chipseq_best_motifs_with_families.tsv
