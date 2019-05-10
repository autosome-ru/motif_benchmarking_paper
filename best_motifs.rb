require 'json'
require_relative 'aucs'
require_relative 'tf_motifs_mapping'

TFCLASS_LEVELS = [:tf_superclass, :tf_class, :tf_family, :tf_subfamily, :tf_genus, :tf_molecular_species]

class Array
  def mean; empty? ? nil : sum(0.0) / size; end
  def median
    return nil if empty?
    sorted = self.sort
    if size % 2 == 1
      sorted[size / 2]
    else
      (sorted[size / 2] + sorted[size / 2 - 1]) / 2.0
    end
  end
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
tfclass_level_name = TFCLASS_LEVELS[tfclass_level - 1]

aucs = Aucs.from_file(aucs_matrix_fn)

header = ['experiment_TF', 'experiment_TF_family', 'best_motif', 'best_auc', 'tfs_of_best_motif', 'best_motifs_family']
puts header.join("\t")
aucs.experiments.group_by{|experiment|
  TF_BY_EXPERIMENT[experiment]
}.each{|experiment_tf, experiments|
  # For each motif we aggregate AUCs over several datasets of an experiment TF
  # We calculate it for all motifs of all TFs, not only an experiment TF
  tf_aucs_by_motif = aucs.motifs.map{|motif|
    motif_aucs = experiments.map{|experiment|
      aucs.auc(motif, experiment)
    }
    [motif, motif_aucs]
  }
  # and select best motif by average AUC over all datasets for a TF
  (best_motif, best_auc) = tf_aucs_by_motif.map{|motif, motif_aucs|
    [motif, motif_aucs.mean]
  }.max_by{|motif, score|
    score
  }

  families_of_experiment = TF_INFO_BY_NAME[experiment_tf][tfclass_level_name]
  tfs_of_best_motif = TFS_BY_MOTIF[best_motif]
  families_of_best_motif = tfs_of_best_motif.flat_map{|tf|
    TF_INFO_BY_NAME[tf][tfclass_level_name]
  }.uniq

  row = [
    experiment_tf, certain_or_ambiguous(families_of_experiment),
    best_motif, best_auc,
    tfs_of_best_motif.join(':'), certain_or_ambiguous(families_of_best_motif),
  ]
  puts row.join("\t")
}
