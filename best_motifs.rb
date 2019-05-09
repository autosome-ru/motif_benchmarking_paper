require 'json'
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
experiment_mapping_fn = ARGV[1]
tfclass_level = Integer(ARGV[2])
tfclass_level_name = TFCLASS_LEVELS[tfclass_level - 1]

motif_tf_pairs = [
  'source_data/motifs/hocomoco_genes2mat.txt',
  'source_data/motifs/jaspar_genes2mat.txt',
].flat_map{|fn|
  File.readlines(fn).map{|l|
    l.split("\t").map(&:strip)
  }
}

motif_tfs = motif_tf_pairs.group_by{|tf, motif|
  motif
}.map{|motif, pairs|
  [motif, pairs.map{|tf, motif| tf }]
}

tf_motifs = motif_tf_pairs.group_by{|tf, motif|
  tf
}.map{|tf, pairs|
  [tf, pairs.map{|tf, motif| motif }]
}

tf_infos = File.readlines('all_tf_infos.json').map{|l|
  JSON.parse(l, symbolize_names: true)
}

TF_INFO_BY_NAME = tf_infos.group_by{|info|
  info[:tf_gene_name]
}.map{|k,vs|
  [k, vs.first]
}.to_h


TFS_BY_MOTIF = motif_tfs.to_h
MOTIFS_BY_TF = tf_motifs.to_h
EXPERIMENT_BY_TF = File.readlines(experiment_mapping_fn).map{|l| l.chomp.split("\t").reverse }.to_h


def tf_by_experiment(experiment)
  EXPERIMENT_BY_TF[experiment]
end

def motif_tf(motif)
  TF_BY_MOTIF[motif]
end

class Aucs
  attr_reader :experiments, :aucs_by_motif
  def initialize(experiments, aucs_by_motif)
    @experiments = experiments
    @experiment_index = @experiments.each_with_index.to_h
    @aucs_by_motif = aucs_by_motif
  end
  def motifs; @aucs_by_motif.keys; end
  def auc(motif, experiment); @aucs_by_motif[motif][ @experiment_index[experiment] ]; end
  def motif_aucs(motif); @experiments.zip(@aucs_by_motif[motif]).to_h; end
  def experiment_aucs(experiment); motifs.map{|motif| [motif, auc(motif, experiment)] }.to_h; end
  def self.from_file(fn)
    File.open(fn){|f|
      experiments = f.readline.chomp.split("\t") # yes, there is no empty cell in the left top corner
      aucs_by_motif = f.readlines.map{|l|
        motif, *aucs = l.chomp.split("\t")
        [motif, aucs.map(&:to_f)]
      }.to_h
      self.new(experiments, aucs_by_motif)
    }
  end
end

aucs = Aucs.from_file(aucs_matrix_fn)
# aucs_unfiltered = Aucs.from_file(aucs_matrix_fn)
# aucs_by_motif_filtered = aucs_unfiltered.aucs_by_motif.reject{|motif,aucs|
#   ["EVX1_HUMAN.H11MO.0.D", "EVX2_HUMAN.H11MO.0.A", "NFKB1_HUMAN.H11MO.0.A"].include?(motif)
# }
# aucs = Aucs.new(aucs_unfiltered.experiments, aucs_by_motif_filtered)

header = ['experiment_TF', 'experiment_TF_family', 'best_motif', 'best_auc', 'tfs_of_best_motif', 'best_motifs_family']
puts header.join("\t")
aucs.experiments.group_by{|experiment|
  tf_by_experiment(experiment)
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
