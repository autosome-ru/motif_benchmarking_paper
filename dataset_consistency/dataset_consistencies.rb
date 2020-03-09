require_relative 'matrix'
require_relative 'distances'
require_relative 'aggregation'
require_relative 'annotations'

##############################

jy_matrix = Matrix.read_matrix('source_data/final/jy10_all_roc.txt', skip_header_cell: false).datasets_renamed{|ds| ds.gsub('-', '.') }.motifs_renamed{|mot| mot.gsub(/^(MA[\d.]+)_.+$/, '\1') }
jolma_matrix = Matrix.read_matrix('source_data/selex/motifs_vs_selex10.tsv', skip_header_cell: false).datasets_renamed{|ds| ds.gsub('-', '.') }.motifs_renamed{|mot| mot.gsub(/^(MA[\d.]+)_.+$/, '\1') }
jy_matrix, jolma_matrix = match_matrices(jy_matrix, jolma_matrix)

##############################

matrix = jolma_matrix

aggregation_level = :tf_family

dataset_groups_none = matrix.datasets.map{|dataset| [dataset, [dataset]] }.sort_by{|ds, _| families_idx_by_dataset(ds, aggregation_level) }.to_h
dataset_groups_tf = dataset_groups_by_tf(matrix.datasets).sort_by{|tf, _| families_idx_by_tf(tf, aggregation_level) }.to_h
dataset_groups_fam = dataset_groups_by_family(matrix.datasets, aggregation_level).sort_by{|fam, _| family_idx(fam) }.to_h

motif_groups_tf = motif_groups_by_tf(matrix.motifs).sort_by{|tf, _| families_idx_by_tf(tf, aggregation_level) }.to_h
motif_groups_fam = motif_groups_by_family(matrix.motifs, aggregation_level).sort_by{|fam, _| family_idx(fam) }.to_h

##############################
mat = aggregate(jolma_matrix, dataset_groups_none, motif_groups_fam)
mat = aggregate(jolma_matrix, dataset_groups_tf, motif_groups_fam)
mat.datasets.each{|dataset|
  # native_fams = mat.motifs & families_by_dataset(dataset, aggregation_level)
  native_fams = mat.motifs & families_by_tf(dataset, aggregation_level)
  foreign_fams = mat.motifs - native_fams
  native_aucs = native_fams.map{|fam| mat.auc(dataset, fam) }
  foreign_aucs = foreign_fams.map{|fam| mat.auc(dataset, fam) }
  if foreign_aucs.max && native_aucs.max && (foreign_aucs.max > native_aucs.max)
    p [dataset, native_fams, native_aucs.max.round(3), foreign_aucs.max.round(3), foreign_fams.max_by{|fam| mat.auc(dataset, fam) || 0 }]
    # puts [native_fams.first, foreign_fams.max_by{|fam| mat.auc(dataset, fam) || 0 }].join("\t")
  else
    # p [dataset, native_fams, native_aucs.max&.round(3)]
  end
}

File.write("dataset_consistency/jolma_none_fam_aggregated.tsv", aggregate(jolma_matrix, dataset_groups_none, motif_groups_fam))
File.write("dataset_consistency/jolma_tf_tf_aggregated.tsv", aggregate(jolma_matrix, dataset_groups_tf, motif_groups_tf))
File.write("dataset_consistency/jolma_fam_tf_aggregated.tsv", aggregate(jolma_matrix, dataset_groups_fam, motif_groups_tf))
File.write("dataset_consistency/jolma_tf_fam_aggregated.tsv", aggregate(jolma_matrix, dataset_groups_tf, motif_groups_fam))
File.write("dataset_consistency/jolma_fam_fam_aggregated.tsv", aggregate(jolma_matrix, dataset_groups_fam, motif_groups_fam))

# tf_jolma_matrix = aggregate(jolma_matrix, dataset_groups_tf, motif_groups_tf)
# tf_jy_matrix = aggregate((jy_matrix, dataset_groups_tf, motif_groups_tf)

# tfs = tf_jolma_matrix.motifs & tf_jolma_matrix.datasets
# deltas = tfs.map{|tf|
#   best_auc_foreign_matrix = (tfs - [tf]).map{|tf_2|
#     tf_jy_matrix.auc(tf,tf_2)
#   }.max
#   delta = best_auc_foreign_matrix - tf_jy_matrix.auc(tf,tf)
#   [tf, delta.round(3)]
# }.sort_by{|k,v| v }.to_h

# motif_consistencies = consistencies(jy_matrix, jolma_matrix)
# tf_consistencies = consistencies(tf_jy_matrix, tf_jolma_matrix)

# ######################

# File.open('/home/ilya/tf_consistencies.tsv', 'w') {|fw|
#   header = ['dataset', 'normed_l1', 'normed_l2', 'l_inf', 'cosine']
#   fw.puts(header.join("\t"))
#   tf_consistencies.each{|ds, info|
#     info = [ds, *info.values_at(:normed_l1, :normed_l2, :l_inf, :cosine)]
#     fw.puts(info.join("\t"))
#   }
# }
