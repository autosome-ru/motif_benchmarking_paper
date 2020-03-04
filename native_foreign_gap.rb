require_relative 'aucs_matrix'
require_relative 'auc_collection'

class Array
  def mean; empty? ? nil : sum(0.0) / size; end
end

def print_matrix_to_file(matrix, col_names:, row_names:, filename:)
  File.open(filename, 'w'){|fw| print_matrix(matrix, col_names: col_names, row_names: row_names, stream: fw) }
end

def print_matrix(matrix, col_names:, row_names:, stream: $stdout)
  stream.puts(['-', *col_names].join("\t"))
  matrix.zip(row_names).each{|row, row_name|
    stream.puts([row_name, *row].join("\t"))
  }
end


def experiment_stats_fam_avemotif(aucs_matrix, annotation, tfclass_level, output_matrix_fn)
  motif_fams = aucs_matrix.motifs.flat_map{|motif| annotation.tfs_by_motif(motif) }.flat_map{|tf| annotation.tf_info_by_gene_name(tf)[tfclass_level] }.uniq #.sort_by{|fam| fam[/\{(.+)\}/,1].split('.').map(&:to_i) }
  experiment_fam_matrix = aucs_matrix.experiments.map{|experiment|
    motif_fams.map{|motif_fam|
      fam_motifs = annotation.motifs_by_tfclass_name(tfclass_level, motif_fam)
      fam_motifs.map{|motif| aucs_matrix.auc(motif, experiment) }.mean
    }
  }

  print_matrix_to_file(experiment_fam_matrix, row_names: aucs_matrix.experiments, col_names: motif_fams, filename: output_matrix_fn)

  experiment_fam_matrix.zip(aucs_matrix.experiments).map{|fam_aucs, experiment|
    auc_by_motif_fam = motif_fams.zip(fam_aucs).to_h
    tf = annotation.tf_by_experiment(experiment)
    native_fams = annotation.tf_info_by_gene_name(tf)[tfclass_level] & motif_fams
    foreign_fams = motif_fams - native_fams

    native_aucs = native_fams.map{|motif_fam| auc_by_motif_fam[motif_fam] }
    foreign_aucs = foreign_fams.map{|motif_fam| auc_by_motif_fam[motif_fam] }

    max_native_auc = native_aucs.max
    max_foreign_auc = foreign_aucs.max

    [experiment, {max_native_auc: max_native_auc, max_foreign_auc: max_foreign_auc, gap: (max_native_auc - max_foreign_auc rescue 0) }]
  }.to_h
end


def experiment_stats_avemotif(aucs_matrix, annotation)
  motif_tfs = aucs_matrix.motifs.flat_map{|motif| annotation.tfs_by_motif(motif) }.uniq
  experiment_tf_matrix = aucs_matrix.experiments.map{|experiment|
    motif_tfs.map{|motif_tf|
      tf_motifs = annotation.motifs_by_tf(motif_tf)
      tf_motifs.map{|motif| aucs_matrix.auc(motif, experiment) }.mean
    }
  }

  experiment_tf_matrix.zip(aucs_matrix.experiments).map{|tf_aucs, experiment|
    auc_by_motif_tf = motif_tfs.zip(tf_aucs).to_h
    tf = annotation.tf_by_experiment(experiment)
    tf_fams = annotation.tf_info_by_gene_name(tf)[:tf_family]
    
    native_tfs = tf_fams.flat_map{|fam| annotation.tfs_by_tfclass_name(:tf_family, fam) }.uniq & motif_tfs
    foreign_tfs = motif_tfs - native_tfs

    native_aucs = native_tfs.map{|motif_tf| auc_by_motif_tf[motif_tf] }
    foreign_aucs = foreign_tfs.map{|motif_tf| auc_by_motif_tf[motif_tf] }

    max_native_auc = native_aucs.max
    max_foreign_auc = foreign_aucs.max

    [experiment, {max_native_auc: max_native_auc, max_foreign_auc: max_foreign_auc, gap: (max_native_auc - max_foreign_auc rescue 0) }]
  }.to_h
end


def experiment_stats(aucs_matrix, annotation, motifs_by_fam)
  aucs_matrix.experiments.map{|experiment|
    tf = annotation.tf_by_experiment(experiment)
    tf_fams = annotation.tf_info_by_gene_name(tf)[:tf_family]

    native_motifs = tf_fams.flat_map{|fam| motifs_by_fam[fam] }.uniq
    native_aucs = native_motifs.map{|motif| aucs_matrix.auc(motif, experiment) }
    max_native_auc = native_aucs.max

    foreign_motifs = aucs_matrix.motifs - native_motifs
    foreign_aucs = foreign_motifs.map{|motif| aucs_matrix.auc(motif, experiment) }
    max_foreign_auc = foreign_aucs.max

    # log_info = [
    #   experiment,
    #   max_native_auc, max_foreign_auc,
    #   native_motifs.map{|motif| [motif, aucs_matrix.auc(motif, experiment)] }.sort_by{|k,v| v }.reverse.first(3).to_h,
    #   foreign_motifs.map{|motif| [motif, aucs_matrix.auc(motif, experiment)] }.select{|motif, v| max_native_auc && v > max_native_auc }.sort_by{|k,v| v }.reverse.first(3).to_h,
    # ]
    # p log_info

    [experiment, {max_native_auc: max_native_auc, max_foreign_auc: max_foreign_auc, gap: (max_native_auc - max_foreign_auc rescue 0) }]
  }.to_h
end

jy_aucs_matrix = AucsMatrix.from_file('source_data/final/jy10_all_roc.txt')
jolma_aucs_matrix = AucsMatrix.from_file('source_data/selex/motifs_vs_selex10.tsv')
annotation = Annotation.new(jy_aucs_matrix.experiments & jolma_aucs_matrix.experiments, jy_aucs_matrix.motifs & jolma_aucs_matrix.motifs)
jy_aucs = AucCollection.new(jy_aucs_matrix, annotation)
jolma_aucs = AucCollection.new(jolma_aucs_matrix, annotation)

families = annotation.tfclass_names(:tf_family)

motifs_by_fam = families.map{|fam|
  [fam, annotation.motifs_by_tfclass_name(:tf_family, fam)]
}.reject{|fam, mots|
  mots.empty?
}.to_h

# jy_stats = experiment_stats(jy_aucs, annotation, motifs_by_fam)
# jolma_stats = experiment_stats(jolma_aucs, annotation, motifs_by_fam)
jy_stats = experiment_stats_fam_avemotif(jy_aucs, annotation, :cisbp_families, '/home/ilya/jy_exp_fam.tsv')
jolma_stats = experiment_stats_fam_avemotif(jolma_aucs, annotation, :cisbp_families, '/home/ilya/jolma_exp_fam.tsv')

experiments = (jy_stats.keys & jolma_stats.keys).sort
# experiments.map{|experiment|
#   [experiment, jolma_stats[experiment][:gap].round(3), jy_stats[experiment][:gap].round(3)]
# }.select{|exp, jolma_gap, jy_gap|
#   jolma_gap < -0.05 || jy_gap < -0.05
# }.sort_by{|exp, jolma_gap, jy_gap| [jy_gap - jolma_gap].min }

# experiments.map{|experiment|
#   [experiment, jolma_stats[experiment].transform_values{|v| v&.round(3) }, jy_stats[experiment].transform_values{|v| v&.round(3) }]
# }.select{|experiment, jolma_info, jy_info|
#   [jolma_info, jy_info].all?{|info| info[:max_native_auc] && info[:max_foreign_auc] }
# }.sort_by{|experiment, jolma_info, jy_info|
#   (jolma_info[:max_foreign_auc] - jy_info[:max_foreign_auc]).abs + (jolma_info[:max_native_auc] - jy_info[:max_native_auc]).abs 
# }

File.open('/home/ilya/experiment_fam_avemotif_aucs.tsv', 'w'){|fw|
  header = ['experiment', 'jolma_max_native', 'jolma_max_foreign', 'jolma_gap', 'jy_max_native', 'jy_max_foreign', 'jy_gap']
  fw.puts(header.join("\t"))
  experiments.each{|experiment|
    info = [
      experiment,
      jolma_stats[experiment].values_at(:max_native_auc, :max_foreign_auc, :gap),
      jy_stats[experiment].values_at(:max_native_auc, :max_foreign_auc, :gap),
    ]
    fw.puts(info.join("\t"))
  }
}
