mkdir -p  results/remap/  results/jy10/  results/jy50/  results/uniprobe/

ruby mean_aucs_matrix.rb source_data/final/remap_all_roc.txt    cisbp_families  >  results/remap/heatmap_cisbpfams_remap.tsv
ruby mean_aucs_matrix.rb source_data/final/jy50_all_roc.txt     cisbp_families  >  results/jy50/heatmap_cisbpfams_jy50.tsv
ruby mean_aucs_matrix.rb source_data/final/jy10_all_roc.txt     cisbp_families  >  results/jy10/heatmap_cisbpfams_jy10.tsv
ruby mean_aucs_matrix.rb source_data/final/uniprobe_all_cor.txt cisbp_families  >  results/uniprobe/heatmap-cor_cisbpfams_uniprobe.tsv

ruby refine_heatmap.rb results/remap/heatmap_cisbpfams_remap.tsv             cisbp_families > results/remap/heatmap_cisbpfams_remap.filtered.tsv
ruby refine_heatmap.rb results/jy50/heatmap_cisbpfams_jy50.tsv               cisbp_families > results/jy50/heatmap_cisbpfams_jy50.filtered.tsv
ruby refine_heatmap.rb results/jy10/heatmap_cisbpfams_jy10.tsv               cisbp_families > results/jy10/heatmap_cisbpfams_jy10.filtered.tsv
ruby refine_heatmap.rb results/uniprobe/heatmap-cor_cisbpfams_uniprobe.tsv   cisbp_families > results/uniprobe/heatmap-cor_cisbpfams_uniprobe.filtered.tsv

python3.6 plot_heatmap.py results/remap/heatmap_cisbpfams_remap.filtered.tsv  results/remap/heatmap_cisbpfams_remap.filtered.svg
python3.6 plot_heatmap.py results/jy50/heatmap_cisbpfams_jy50.filtered.tsv    results/jy50/heatmap_cisbpfams_jy50.filtered.svg
python3.6 plot_heatmap.py results/jy10/heatmap_cisbpfams_jy10.filtered.tsv    results/jy10/heatmap_cisbpfams_jy10.filtered.svg
python3.6 plot_heatmap.py results/uniprobe/heatmap-cor_cisbpfams_uniprobe.filtered.tsv  results/uniprobe/heatmap-cor_cisbpfams_uniprobe.filtered.svg --correlation
