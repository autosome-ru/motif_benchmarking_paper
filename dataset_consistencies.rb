class Matrix
  attr_reader :motifs, :datasets, :motif_dataset_matrix
  def initialize(datasets, motifs, motif_dataset_matrix)
    @datasets = datasets
    @motifs = motifs
    @motif_dataset_matrix = motif_dataset_matrix
    @dataset_motif_matrix = motif_dataset_matrix.transpose
  end
  def self.read_matrix(filename)
    lns = File.readlines(filename).map{|l| l.chomp.split("\t") }
    datasets = lns.first
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
end

def norm(a)
  a.map{|ai| ai ** 2 }.sum(0.0) ** 0.5
end

def cosine(a, b)
  a.zip(b).map{|ai, bi| ai * bi }.sum(0.0) / (norm(a) * norm(b))
end

def normed_l1(a, b)
  a.zip(b).map{|ai, bi| (ai - bi).abs }.sum(0.0) / a.size
end

def normed_l2(a, b)
  (a.zip(b).map{|ai, bi| (ai - bi) ** 2 }.sum(0.0) / a.size) ** 0.5 
end


def l_inf(a, b)
  a.zip(b).map{|ai, bi| (ai - bi).abs }.max
end


def consistencies(matrix_1, matrix_2)
  datasets = (matrix_1.datasets & matrix_2.datasets)
  motifs = (matrix_1.motifs & matrix_2.motifs)

  matrix_2 = matrix_2.bound_at(datasets, motifs)
  matrix_1 = matrix_1.bound_at(datasets, motifs)

  datasets.map{|dataset|
    jy_vect = matrix_1.dataset_aucs(dataset)
    jolma_vect = matrix_2.dataset_aucs(dataset)
    info = {
      normed_l1: normed_l1(jy_vect, jolma_vect),
      normed_l2: normed_l2(jy_vect, jolma_vect),
      l_inf: l_inf(jy_vect, jolma_vect),
      cosine: cosine(jy_vect, jolma_vect)
    }.map{|k, v|
      [k, v.round(3)]
    }.to_h
    [dataset, info]
  }.to_h
end


def match_matrices(matrix_1, matrix_2)
  datasets = (matrix_1.datasets & matrix_2.datasets)
  motifs = (matrix_1.motifs & matrix_2.motifs)
  [matrix_1, matrix_2].map{|matrix| matrix.bound_at(datasets, motifs) }
end

############################

def tf_motif_pairs
  @tf_motif_pairs ||= [
    ['source_data/motifs/hocomoco_genes2mat.txt', ->(motif_name){ motif_name }],
    ['source_data/motifs/jaspar_genes2mat.txt', ->(motif_name){ motif_name.split('_').first }],
    ['ccg.epfl.ch/pwmtools/benchmarking/genes2cisbp.txt', ->(motif_name){ motif_name }],
  ].flat_map{|fn, motif_name_transformation|
    File.readlines(fn).map{|l|
      l.split("\t").map(&:strip)
    }.map{|tf, motif_name, *rest|
      [tf, motif_name_transformation.call(motif_name)]
    }
  }
end


def tfs_by_motif(motif)
  @tfs_by_motif ||= tf_motif_pairs.group_by{|tf, motif|
    motif
  }.map{|motif, pairs|
    [motif, pairs.map{|tf, motif| tf }]
  }.to_h
  @tfs_by_motif[motif] || []
end

def motifs_by_tf(tf)
  @motifs_by_tf ||= tf_motif_pairs.group_by{|tf, motif|
    tf
  }.map{|tf, pairs|
    [tf, pairs.map{|tf, motif| motif }]
  }.to_h
  @motifs_by_tf[tf] || []
end

############################

def experiment_tf_pairs
  @experiment_tf_pairs ||= [
    'source_data/chipseq/remap_genes2exp.txt',
    'source_data/selex/jolma13_genes2exp.txt',
    'ccg.epfl.ch/pwmtools/benchmarking/genes2jolma_yang.txt',
    'ccg.epfl.ch/pwmtools/benchmarking/genes2uniprobe.txt',
    'source_data/uniprobe/genes2uniprobe_manual.txt',
  ].flat_map{|experiment_mapping_fn|
    File.readlines(experiment_mapping_fn).map{|l|
      tf, experiment = l.chomp.split("\t").first(2)
      [experiment, tf]
    }
  }
end

def tf_by_experiment(experiment)
  @tf_by_experiment ||= experiment_tf_pairs.to_h
  @tf_by_experiment[experiment]
end

def experiments_by_tf(tf)
  @experiments_by_tf ||= experiment_tf_pairs.group_by{|experiment, tf|
    tf
  }.map{|tf, exp_tf_pairs|
    [tf, exp_tf_pairs.map{|exp,tf| exp }]
  }.to_h
  @experiments_by_tf[tf] || []
end

############################

class Array
  def mean; sum(0.0) / length; end
end


def aggregate_matrix(matrix, only_common: false)
  motif_tfs = matrix.motifs.flat_map{|motif| tfs_by_motif(motif) }.uniq
  dataset_tfs = matrix.datasets.map{|dataset| tf_by_experiment(dataset) }.uniq
  common_tfs = motif_tfs & dataset_tfs
  if only_common
    motif_tfs = common_tfs
    dataset_tfs = common_tfs
  end

  new_matrix = motif_tfs.map{|motif_tf|
    tf_motifs = motifs_by_tf(motif_tf).select{|motif| matrix.motifs.include?(motif) }
    dataset_tfs.map{|dataset_tf|
      tf_datasets = experiments_by_tf(dataset_tf).select{|dataset| matrix.datasets.include?(dataset) }
      aucs = tf_motifs.flat_map{|motif|
        tf_datasets.map{|dataset|
          matrix.auc(dataset, motif)
        }
      }
      aucs.mean
    }
  }
  Matrix.new(dataset_tfs, motif_tfs, new_matrix)
end

def aggregate_motifs(matrix)
  motif_tfs = matrix.motifs.flat_map{|motif| tfs_by_motif(motif) }.uniq
  new_matrix = motif_tfs.map{|motif_tf|
    tf_motifs = motifs_by_tf(motif_tf).select{|motif| matrix.motifs.include?(motif) }
    matrix.datasets.map{|dataset|
      aucs = tf_motifs.flat_map{|motif|
        matrix.auc(dataset, motif)
      }
      aucs.mean
    }
  }
  Matrix.new(matrix.datasets, motif_tfs, new_matrix)
end

def aggregate_datasets(matrix)
  dataset_tfs = matrix.datasets.map{|dataset| tf_by_experiment(dataset) }.uniq
  new_matrix = matrix.motifs.map{|motif|
    dataset_tfs.map{|dataset_tf|
      tf_datasets = experiments_by_tf(dataset_tf).select{|dataset| matrix.datasets.include?(dataset) }
      aucs = tf_datasets.map{|dataset|
        matrix.auc(dataset, motif)
      }
      aucs.mean
    }
  }
  Matrix.new(dataset_tfs, matrix.motifs, new_matrix)
end

def leave_common(matrix)
  dataset_tfs = matrix.datasets.map{|dataset| tf_by_experiment(dataset) }.uniq
  motif_tfs = matrix.motifs.flat_map{|motif| tfs_by_motif(motif) }.uniq
  common_tfs = motif_tfs & dataset_tfs
  datasets = common_tfs.flat_map{|tf| experiments_by_tf(tf) }.select{|dataset| matrix.datasets.include?(dataset) }.uniq
  motifs = common_tfs.flat_map{|tf| motifs_by_tf(tf) }.select{|motif| matrix.motifs.include?(motif) }.uniq
  matrix.bound_at(datasets, motifs)
end

##############################


jy_matrix = Matrix.read_matrix('source_data/final/jy10_all_roc.txt').datasets_renamed{|ds| ds.gsub('-', '.') }.motifs_renamed{|mot| mot.gsub(/^(MA[\d.]+)_.+$/, '\1') }
jolma_matrix = Matrix.read_matrix('source_data/selex/motifs_vs_selex10.tsv').datasets_renamed{|ds| ds.gsub('-', '.') }.motifs_renamed{|mot| mot.gsub(/^(MA[\d.]+)_.+$/, '\1') }
jy_matrix, jolma_matrix = match_matrices(jy_matrix, jolma_matrix)

tf_jolma_matrix = aggregate_datasets(aggregate_motifs(jolma_matrix))
tf_jy_matrix = aggregate_datasets(aggregate_motifs(jy_matrix))

tfs = tf_jolma_matrix.motifs & tf_jolma_matrix.datasets
deltas = tfs.map{|tf|
  best_auc_foreign_matrix = (tfs - [tf]).map{|tf_2|
    tf_jy_matrix.auc(tf,tf_2)
  }.max
  delta = best_auc_foreign_matrix - tf_jy_matrix.auc(tf,tf)
  [tf, delta.round(3)]
}.sort_by{|k,v| v }.to_h



motif_consistencies = consistencies(jy_matrix, jolma_matrix)
tf_consistencies = consistencies(tf_jy_matrix, tf_jolma_matrix)

######################

File.open('/home/ilya/tf_consistencies.tsv', 'w') {|fw|
  header = ['dataset', 'normed_l1', 'normed_l2', 'l_inf', 'cosine']
  fw.puts(header.join("\t"))
  tf_consistencies.each{|ds, info|
    info = [ds, *info.values_at(:normed_l1, :normed_l2, :l_inf, :cosine)]
    fw.puts(info.join("\t"))
  }
}



