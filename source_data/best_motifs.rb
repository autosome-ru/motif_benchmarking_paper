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

TF_BY_MOTIF = File.readlines('annotation/motif_annotation.tsv').map{|l| l.chomp.split("\t", 4).first(2) }.to_h
FAMILY_BY_TF = File.readlines('annotation/tf_families.tsv').map{|l| l.chomp.split("\t", 3).first(2) }.to_h
EXPERIMENT_BY_TF = File.readlines('remap_genes2exp.txt').map{|l| l.chomp.split("\t").reverse }.to_h
GENE_BY_UNIPROT = File.readlines('annotation/hoco_ann.tsv').map{|l|
  motif, gene = l.chomp.split("\t").first(2)
  [motif.split('.').first, gene]
}.to_h

def tf_by_experiment(experiment)
  # experiment.split('.')[1]
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
      experiments = f.readline.chomp.split("\t")
      aucs_by_motif = f.readlines.map{|l|
        motif, *aucs = l.chomp.split("\t")
        [motif, aucs.map(&:to_f)]
      }.to_h
      self.new(experiments, aucs_by_motif)
    }
  end
end

# hocomoco(chipseq)
aucs_unfiltered = Aucs.from_file('motifs_vs_remap.tsv')
aucs_by_motif_filtered = aucs_unfiltered.aucs_by_motif.reject{|motif,aucs|
  ["EVX1_HUMAN.H11MO.0.D", "EVX2_HUMAN.H11MO.0.A", "NFKB1_HUMAN.H11MO.0.A"].include?(motif)
}
# pp aucs_unfiltered.motifs
aucs = Aucs.new(aucs_unfiltered.experiments, aucs_by_motif_filtered)
# aucs.experiments.first(50).map do |experiment|

#   # best_motifs = aucs.experiment_aucs(experiment).sort_by{|motif, auc|
#   #   auc
#   # }.reverse.first(5).to_h

#   tf_aucs = aucs.experiment_aucs(experiment).map{|motif, auc|
#     [motif_tf(motif), auc]
#   }.group_by{|tf,auc| tf }.map{|tf, auc_group|
#     [tf, auc_group.map{|tf, auc| auc }]
#   }
#   infos = tf_aucs.map{|tf, aucs| [GENE_BY_UNIPROT[tf], aucs.max] }.sort_by(&:last).reverse
  
#   self_tf = tf_by_experiment(experiment)
#   _, self_auc = infos.detect{|tf, auc| tf == self_tf }
#   best_tf, best_auc = infos.first
#   # puts [experiment, best_tf, best_auc, self_auc, self_auc && (best_auc - self_auc)].join("\t")
# end

motifs_by_tf = aucs.motifs.group_by{|motif| motif_tf(motif) }.to_h


# header = [
#   'TF',
#   [['best_mean_motif', 'best_mean_auc'], ['best_proper_mean', 'best_proper_mean_auc'], 'mean_delta'],
#   [['best_median_motif', 'best_median_auc'], ['best_proper_median', 'best_proper_median_auc'], 'median_delta'],
#   [['best_max_motif', 'best_max_auc'], ['best_proper_max', 'best_proper_max_auc'], 'max_delta'],
# ]
# puts header.join("\t")
# aucs.experiments.group_by{|experiment|
#   tf_by_experiment(experiment)
# }.each{|experiment_tf, experiments|
#   # for datasets of certain TF we test motifs of all factors
#   tf_aucs_by_motif = aucs.motifs.map{|motif|
#     motif_aucs = experiments.map{|experiment|
#       aucs.auc(motif, experiment)
#     }
#     [motif, motif_aucs]
#   }
#   # and select best motif by average AUC over all datasets for a TF
#   best_mean = tf_aucs_by_motif.map{|motif, motif_aucs| [motif, motif_aucs.mean] }.max_by{|motif, score| score }
#   best_median = tf_aucs_by_motif.map{|motif, motif_aucs| [motif, motif_aucs.median] }.max_by{|motif, score| score }
#   best_max = tf_aucs_by_motif.map{|motif, motif_aucs| [motif, motif_aucs.max] }.max_by{|motif, score| score }

#   proper_motifs = motifs_by_tf[experiment_tf]
#   if proper_motifs
#     best_proper_mean = tf_aucs_by_motif.select{|motif, _| proper_motifs.include?(motif) }.map{|motif, motif_aucs| [motif, motif_aucs.mean] }.max_by(&:last)
#     best_proper_median = tf_aucs_by_motif.select{|motif, _| proper_motifs.include?(motif) }.map{|motif, motif_aucs| [motif, motif_aucs.median] }.max_by(&:last)
#     best_proper_max = tf_aucs_by_motif.select{|motif, _| proper_motifs.include?(motif) }.map{|motif, motif_aucs| [motif, motif_aucs.max] }.max_by(&:last)

#     mean_delta = best_mean.last - best_proper_mean.last
#     median_delta = best_median.last - best_proper_median.last
#     max_delta = best_max.last - best_proper_max.last
#     puts [experiment_tf, [best_mean, best_proper_mean, mean_delta], [best_median, best_proper_median, median_delta], [best_max, best_proper_max, max_delta]].join("\t")
#   else
#     puts [experiment_tf, [best_mean, [nil, nil], nil], [best_median, [nil,nil], nil], [best_max, [nil, nil], nil]].join("\t")
#   end
# }


# header = ['family', 'best_mean_motif', 'best_mean_auc',]
# puts header.join("\t")
# aucs.experiments.group_by{|experiment|
#   FAMILY_BY_TF[tf_by_experiment(experiment)] || ''
# }.each{|family, experiments|
#   # for datasets of certain TF we test motifs of all factors
#   tf_aucs_by_motif = aucs.motifs.map{|motif|
#     motif_aucs = experiments.map{|experiment|
#       aucs.auc(motif, experiment)
#     }
#     [motif, motif_aucs]
#   }
#   # and select best motif by average AUC over all datasets for a TF
#   best_mean = tf_aucs_by_motif.map{|motif, motif_aucs| [motif, motif_aucs.mean] }.max_by{|motif, score| score }
#   best_mean_motif, best_mean_auc = *best_mean
#   best_mean_tf = TF_BY_MOTIF[best_mean_motif]

#   puts [family, best_mean_motif, best_mean_auc, best_mean_tf].join("\t")
# }


header = ['TF', 'best_mean_motif', 'best_mean_auc', 'best_mean_tf']
puts header.join("\t")
aucs.experiments.group_by{|experiment|
  tf_by_experiment(experiment)
}.each{|experiment_tf, experiments|
  # for datasets of certain TF we test motifs of all factors
  tf_aucs_by_motif = aucs.motifs.map{|motif|
    motif_aucs = experiments.map{|experiment|
      aucs.auc(motif, experiment)
    }
    [motif, motif_aucs]
  }
  # and select best motif by average AUC over all datasets for a TF
  best_mean = tf_aucs_by_motif.map{|motif, motif_aucs| [motif, motif_aucs.mean] }.max_by{|motif, score| score }
  best_mean_motif, best_mean_auc = *best_mean
  best_mean_tf = TF_BY_MOTIF[best_mean_motif]
  puts [experiment_tf, best_mean_motif, best_mean_auc, best_mean_tf].join("\t")
}



# TFInfo = Struct.new(:tf, :exp, :num_ds, :winning_motif, :winning_tf, :auc) do
#   def to_s; [tf, exp, num_ds, winning_motif, winning_tf, auc].join("\t"); end
# end
# best_motifs = File.readlines('all_best_for_tf.txt').map{|l|
#   tf, exp, num_ds, winning_motif, winning_tf, auc = l.chomp.split("\t")
#   TFInfo.new(tf, exp, num_ds.to_i, winning_motif, winning_tf, auc.to_f)
# }
