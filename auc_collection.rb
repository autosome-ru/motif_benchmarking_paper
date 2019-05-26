require_relative 'aucs_matrix'
require_relative 'annotation'

class AucCollection
  attr_reader :annotation, :aucs_matrix
  def initialize(aucs_matrix, annotation)
    @aucs_matrix = aucs_matrix # AucsMatrix.from_file('source_data/chipseq/motifs_vs_remap.tsv')
    @annotation = annotation
  end

  def auc(motif, experiment); @aucs_matrix.auc(motif, experiment); end
  def motifs; @aucs_matrix.motifs; end
  def experiments; @aucs_matrix.experiments; end

  def aucs_subset(motifs, experiments)
    motifs = motifs.select{|motif| @aucs_matrix.motifs.include?(motif) }
    experiments = experiments.select{|experiment| @aucs_matrix.experiments.include?(experiment) }
    motifs.flat_map{|motif|
      experiments.map{|experiment|
        aucs_matrix.auc(motif, experiment)
      }
    }
  end

  # def aucs_over_tf_datasets(motifs, experiment_tfs)
  #   experiment_datasets = experiment_tfs.flat_map{|tf|
  #     experiments_by_tf(tf)
  #   }.uniq
  #   aucs_subset(motifs, experiment_datasets)
  # end

  def aucs_over_tfclass_group_datasets(motifs, tfclass_level, tfclass_name)
    experiment_datasets = annotation.experiments_by_tfclass_name(tfclass_level, tfclass_name).uniq
    aucs_subset(motifs, experiment_datasets)
  end
end
