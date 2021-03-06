require 'json'
require_relative 'aucs_matrix'
require_relative 'annotation'

class Array
  def mean; empty? ? nil : sum(0.0) / size; end
end

def certain_or_ambiguous(values)
  if values.size == 1
    values.first
  elsif values.size > 1
    'ambiguous'
  elsif values.size == 0
    'unknown'
  end
end

only_common_tf_experiments = ARGV.delete('--only-common-experiment-tfs')

aucs_matrix_fn = ARGV[0]
tfclass_level_name = ARGV[1].to_sym
auc_threshold = Float(ARGV[2])

aucs_matrix = AucsMatrix.from_file(aucs_matrix_fn)
annotation = Annotation.new(aucs_matrix.experiments, aucs_matrix.motifs)

experiment_tfs = annotation.experiment_tfs

if only_common_tf_experiments
  aucs_matrix_chipseq = AucsMatrix.from_file('source_data/final/remap_all_roc.txt')
  aucs_matrix_selex = AucsMatrix.from_file('source_data/final/jy10_all_roc.txt')
  chipseq_annotation = Annotation.new(aucs_matrix_chipseq.experiments, aucs_matrix_chipseq.motifs)
  selex_annotation = Annotation.new(aucs_matrix_selex.experiments, aucs_matrix_selex.motifs)

  common_experiment_tfs = (chipseq_annotation.experiment_tfs + selex_annotation.experiment_tfs).uniq.select{|tf|
    chipseq_annotation.experiments_by_tf(tf).size >= 1 &&
    selex_annotation.experiments_by_tf(tf).size >= 1
  }
  experiment_tfs = experiment_tfs & common_experiment_tfs
end

header = ['experiment_TF', 'experiment_family', 'num_experiments', 'motif', 'auc', 'tfs_of_motif', 'motif_family']
puts header.join("\t")
experiment_tfs.flat_map{|experiment_tf|
  experiments = annotation.experiments_by_tf(experiment_tf)
  # For each motif we aggregate AUCs over several datasets of an experiment TF
  # We calculate it for all motifs of all TFs, not only an experiment TF
  tf_aucs_by_motif = aucs_matrix.motifs.map{|motif|
    motif_aucs = experiments.map{|experiment|
      aucs_matrix.auc(motif, experiment)
    }
    [motif, motif_aucs]
  }
  tf_aucs_by_motif.map{|motif, aucs|
    {experiment_tf: experiment_tf, num_experiments: experiments.size, motif: motif, auc: aucs.mean}
  }
}.select{|info|
  info[:auc] >= auc_threshold;
}.each{|info|
  tf_info = annotation.tf_info_by_gene_name(info[:experiment_tf])
  next  unless tf_info
  families_of_experiment = tf_info[tfclass_level_name]
  tfs_of_motif = annotation.tfs_by_motif(info[:motif])
  families_of_motif = tfs_of_motif.flat_map{|tf|
    annotation.tf_info_by_gene_name(tf)[tfclass_level_name]
  }.uniq

  row = [
    info[:experiment_tf], certain_or_ambiguous(families_of_experiment), info[:num_experiments],
    info[:motif], info[:auc],
    tfs_of_motif.join(':'), certain_or_ambiguous(families_of_motif),
  ]
  puts row.join("\t")
}
