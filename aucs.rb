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
  def motif_aucs_across_experiments(motif, experiments_subset)
    experiments_subset.map{|experiment|
      [experiment, auc(motif, experiment)]
    }.to_h
  end
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
