require 'json'

class Matrix
  attr_reader :motifs, :datasets, :motif_dataset_matrix
  def initialize(datasets, motifs, motif_dataset_matrix)
    @datasets = datasets
    @motifs = motifs
    @motif_dataset_matrix = motif_dataset_matrix
    @dataset_motif_matrix = motif_dataset_matrix.transpose
  end

  def self.read_matrix(filename, skip_header_cell: true)
    lns = File.readlines(filename).map{|l| l.chomp.split("\t") }
    datasets = lns.first
    datasets = datasets.drop(1)  if skip_header_cell
    motifs = lns.drop(1).map(&:first)
    motif_dataset_matrix = lns.drop(1).map{|l| l.drop(1).map(&:to_f) }
    self.new(datasets, motifs, motif_dataset_matrix)
  end

  def motif_index(motif)
    @motif_index_cache ||= @motifs.each_with_index.to_h
    @motif_index_cache[motif]
  end

  def dataset_index(dataset)
    @dataset_index_cache ||= @datasets.each_with_index.to_h
    @dataset_index_cache[dataset]
  end

  def motif_aucs(motif)
    @motif_dataset_matrix[ motif_index(motif) ]
  end

  def dataset_aucs(dataset)
    @dataset_motif_matrix[ dataset_index(dataset) ]
  end

  def auc(dataset, motif)
    @motif_dataset_matrix[ motif_index(motif) ][ dataset_index(dataset) ]
  end

  def bound_at(datasets_subset, motifs_subset)
    new_matrix = motifs_subset.map{|motif|
      datasets_subset.map{|dataset|
        auc(dataset, motif)
      }
    }
    self.class.new(datasets_subset, motifs_subset, new_matrix)
  end

  def motifs_renamed(&block)
    self.class.new(datasets, motifs.map(&block), motif_dataset_matrix)
  end

  def datasets_renamed(&block)
    self.class.new(datasets.map(&block), motifs, motif_dataset_matrix)
  end

  def to_s(skip_header_cell: true)
    header = skip_header_cell ? ['-', *datasets] : datasets
    [
      header,
      *motifs.zip(motif_dataset_matrix).map{|motif, motif_row|
        [motif, *motif_row]
      }
    ].map{|row| row.join("\t") }.join("\n")
  end
end

def match_matrices(matrix_1, matrix_2)
  datasets = (matrix_1.datasets & matrix_2.datasets)
  motifs = (matrix_1.motifs & matrix_2.motifs)
  [matrix_1, matrix_2].map{|matrix| matrix.bound_at(datasets, motifs) }
end
