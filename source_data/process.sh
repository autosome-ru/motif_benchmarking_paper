cat jaspar_genes2mat.txt | tr -d '\r' | sed -re 's/ +/\t/g' | sponge jaspar_genes2mat.txt
cat hocomoco_genes2mat.txt | tr -d '\r' | sed -re 's/ +/\t/g' | sponge hocomoco_genes2mat.txt
cat remap_genes2exp.txt | tr -d '\r' | sed -re 's/ +/\t/g' | sponge remap_genes2exp.txt

cat hocomoco_remap_all_vs_all.txt | sed -re 's/^ +//' | sed -re 's/ +/\t/g' | sponge hocomoco_remap_all_vs_all.txt
cat jaspar_remap_all_vs_all.txt | sed -re 's/^ +//' | sed -re 's/ +/\t/g' | sponge jaspar_remap_all_vs_all.txt

cut -f1 hocomoco_remap_all_vs_all.txt | tail -n+2 | sort > hocomoco_motifs.txt
cut -f1 jaspar_remap_all_vs_all.txt | tail -n+2 | sort > jaspar_motifs.txt


cp ~/iogen/projects/hocomoco11/final_bundle/hocomoco11/full/HUMAN/mono/pwm/*.pwm all_hoco_motifs
cp ~/iogen/projects/hocomoco11/final_bundle/hocomoco11/full/MOUSE/mono/pwm/*.pwm all_hoco_motifs

mkdir cluster_centers
cut -f1 annotation/hocomoco_cluster_list_192.txt | xargs -n1 -I{} cp all_hoco_motifs/{}.pwm cluster_centers

mkdir core_motifs
cut -f1 annotation/hocomoco_cluster_list_192.txt 
cat annotation/hocomoco_cluster_list_192.txt | tr '\t;' '\n\n' | sort | uniq | xargs -n1 -I{} cp all_hoco_motifs/{}.pwm core_motifs

cat annotation/HOCOMOCOv11_full_annotation_* | fgrep -f hocomoco_motifs.txt | cut -f 1,2,14,15 -d $'\t' | ruby cluster_centers_by_motif.rb | sort -k1,1 > annotation/hoco_ann.tsv
cat annotation/hoco_ann.tsv | awk -F $'\t' -e '($5=="-"){print $1}' | xargs -n1 -I{} echo './nearest_cluster.sh {} cluster_centers' | bash > annotation/missing_cluster_centers.tsv
cat annotation/hoco_ann.tsv | awk -F $'\t' -e '($5=="-"){print $1}' | xargs -n1 -I{} echo './nearest_cluster.sh {} core_motifs' | bash > annotation/missing_cluster_representatives.tsv

ruby complete_annotation_w_clusters.rb > annotation/hoco_ann_refined.tsv 


join jaspar_genes2mat.txt hocomoco_genes2mat.txt -a 1 -t $'\t' | sort -k1,1 > jaspar2hocomoco.tsv
join -a 1 -t $'\t' -1 3 <( sort -k3,3 jaspar2hocomoco.tsv ) annotation/hoco_ann.tsv | cut -d $'\t' -f3- | sort | uniq | cut -d $'\t' -f1-4 > annotation/jasp_ann.tsv

# этим руками придется проверить аннотацию семейств:
cat jaspar2hocomoco.tsv  | grep -v '.H11MO.' | awk -F $'\t' -e '{print $2"\t"$1 }' | sort | uniq | xargs -n1 -I{} echo "echo -n {} $'\t'; ./jaspar_families.sh '{}'" | bash >> annotation/jasp_ann.tsv

# руками удалены мотивы, у которых есть две строки: с семейством и без

# В хокомоке это ретрактнуто: NFKB1_HUMAN, EVX1, EVX2 и я их не возвращал
# В джаспаре это добавлено руками
# MA0109.1_HLTF	HLTF		
# MA0111.1_Spz1	SPZ1		
# MA0637.1_CENPB	CENPB		
# MA0105.4_NFKB1
# MA0887.1_EVX1
# MA0888.1_EVX2

# уберем дубликаты (появились оттого, что одному джаспарному мотиву соответствует несколько хокомочных - и с каждого из них берется семейство, одно и то же)
cat annotation/jasp_ann.tsv | sort | uniq | sort -k1,1 | sponge annotation/jasp_ann.tsv

cat annotation/jasp_ann.tsv | sort | uniq | grep -vPe '^[\w.]+_(NFKB1|EVX1|EVX2)$' | sponge annotation/jasp_ann.tsv

# В hoco_ann_refined мотив MGAP_HUMAN.H11MO.0.D раздвоен, чтобы дать по одной записи на каждое семейство
cat annotation/jasp_ann.tsv | ruby -e 'readlines.map{|l| l.chomp.split("\t") }.group_by{|r| r.first }.map{|motif, rows| genes = rows.map{|r| r[1] }; fams = rows.map{|r| r[2..3] }; [rows.first.first, genes, fams] }.select{|motif, genes, fams| fams.uniq.size == 1 }.each{|motif, genes, fams| puts [motif, genes.uniq.join(";"), fams.uniq].flatten.join("\t") }' > annotation/jasp_ann_single_family.tsv
cat annotation/jasp_ann.tsv | ruby -e 'readlines.map{|l| l.chomp.split("\t") }.group_by{|r| r.first }.map{|motif, rows| genes = rows.map{|r| r[1] }; fams = rows.map{|r| r[2..3] }; [rows.first.first, genes, fams] }.select{|motif, genes, fams| fams.uniq.size != 1 }.each{|motif, genes, fams| genes.zip(fams).each{|gene,fam| puts [motif, gene, fam].flatten.join("\t") }}' > annotation/jasp_ann_several_families.tsv

cat annotation/hoco_ann_refined.tsv annotation/jasp_ann.tsv | cut -d $'\t' -f2,3,4 | sort | uniq | sed -re 's/^(\S+)$/\1\t\t/' > annotation/tf_families.tsv
cat annotation/hoco_ann_refined.tsv annotation/jasp_ann.tsv | cut -d $'\t' -f1,2,3,4 | sort | uniq | sed -re 's/^(\S+)$/\1\t\t/' > annotation/motif_annotation.tsv

cat hocomoco_remap_all_vs_all.txt <( tail -n+2 jaspar_remap_all_vs_all.txt ) > motifs_vs_remap.tsv


join -t $'\t' <( sort -k1,1 remap_genes2exp.txt ) <( sort -k1,1 annotation/tf_families.tsv ) | sort | uniq > annotation/experiment_family.tsv

# здесь появляются лишние пробелы после запроса, их я удалял вручную
cat unknown_families.txt | xargs -n1 -I{} echo "echo -n {} \$'\\t'; curl -qs 'https://www.uniprot.org/uniprot/?query={}+organism%3AHuman&sort=score&format=tab&columns=id,entry%20name,genes&limit=1' | tail -1" | bash > unknown_families_augmented.txt
cat unknown_families_augmented.txt | ruby assign_families.rb 3 | cut -f 1,3,7,8 > unknown_families_recognized.tsv

ruby augment_jaspar_annotation.rb > annotation/jasp_ann_augmented.tsv


# здесь появляются лишние пробелы после запроса, их я удалял вручную
cat annotation/jasp_ann_augmented.tsv | cut -f2 | sort | uniq | xargs -n1 -I{} echo "echo -n {} \$'\\t'; curl -qs 'https://www.uniprot.org/uniprot/?query=gene_exact:{}+organism:Human&sort=score&format=tab&columns=id,entry%20name,genes&limit=1' | tail -1" | bash | cut -f1,3 > annotation/jaspar_uniprots.tsv

join -t $'\t' -2 2 <( sort -k1,1 -t $'\t' annotation/jaspar_uniprots.tsv ) <( sort -t $'\t' -k2,2 annotation/jasp_ann_augmented.tsv ) | awk -F $'\t' -e '{print $3 "\t" $2 "\t" $1 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' > annotation/jaspar_prefinal.tsv

cat annotation/jaspar_prefinal.tsv | ruby assign_families.rb 2 > annotation/jaspar_prefinal_2.tsv
cat annotation/jaspar_prefinal_2.tsv | cut -d $'\t' -f 1,2,3,8,9,10-13 > annotation/jaspar_prefinal_3.tsv 
cat annotation/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv | tail -n+2 | ruby assign_families.rb 18 > annotation/hocomoco_prefinal.tsv
cat annotation/hocomoco_prefinal.tsv | awk -F $'\t' -e '{print $1 "\t" $18 "\t" $2 "\t" $8 "\t" $3 "\t" $20 "\t" $21 "\t" $22 "\t" $23}' > annotation/hocomoco_prefinal_2.tsv 

ruby best_motifs.rb > best_motifs.tsv
cat annotation/jaspar_prefinal_3.tsv annotation/hocomoco_prefinal_2.tsv > annotation/motifs_prefinal.tsv
