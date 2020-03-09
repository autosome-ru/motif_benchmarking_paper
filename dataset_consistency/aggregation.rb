
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

def aggregate_motifs(matrix, motif_groups)
  motif_group_names = motif_groups.keys
  new_matrix = motif_group_names.map{|group|
    motifs = motif_groups[group]
    matrix.datasets.map{|dataset|
      aucs = motifs.map{|motif| matrix.auc(dataset, motif) }
      aucs.mean
    }
  }
  Matrix.new(matrix.datasets, motif_group_names, new_matrix)
end

def aggregate_datasets(matrix, dataset_groups)
  dataset_group_names = dataset_groups.keys
  new_matrix = matrix.motifs.map{|motif|
    dataset_group_names.map{|group|
      datasets = dataset_groups[group]
      aucs = datasets.map{|dataset| matrix.auc(dataset, motif) }
      aucs.mean
    }
  }
  Matrix.new(dataset_group_names, matrix.motifs, new_matrix)
end

def aggregate(matrix, dataset_groups, motif_groups)
  aggregate_datasets(aggregate_motifs(matrix, motif_groups), dataset_groups)
end
