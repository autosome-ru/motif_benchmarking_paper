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

aucs_matrix_fn = ARGV[0]
tfclass_level = Integer(ARGV[1])
auc_threshold = Float(ARGV[2])
tfclass_level_name = TFCLASS_LEVELS[tfclass_level - 1]

aucs_matrix = AucsMatrix.from_file(aucs_matrix_fn)
annotation = Annotation.new(aucs_matrix.experiments, aucs_matrix.motifs)

header = ['experiment_TF', 'experiment_family', 'num_experiments', 'motif', 'auc', 'tfs_of_motif', 'motif_family']
puts header.join("\t")
aucs_matrix.experiments.group_by{|experiment|
  annotation.tf_by_experiment(experiment)
}.flat_map{|experiment_tf, experiments|
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
  info[:auc] >= auc_threshold
}.each{|info|
  families_of_experiment = annotation.tf_info_by_gene_name(info[:experiment_tf])[tfclass_level_name]
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
