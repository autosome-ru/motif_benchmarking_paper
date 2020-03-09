require 'json'

def tfclass_infos
  @tfclass_infos ||= File.readlines('all_tf_infos.json').map{|l|
    JSON.parse(l, symbolize_names: true)
  }
end

def tf_info_by_gene_name(name)
  @tf_info_by_gene_name ||= tfclass_infos.group_by{|info|
    info[:tf_gene_name]
  }.map{|k,vs|
    [k, vs.first]
  }.to_h
  @tf_info_by_gene_name[name]
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

def leave_common(matrix)
  dataset_tfs = matrix.datasets.map{|dataset| tf_by_experiment(dataset) }.uniq
  motif_tfs = matrix.motifs.flat_map{|motif| tfs_by_motif(motif) }.uniq
  common_tfs = motif_tfs & dataset_tfs
  datasets = common_tfs.flat_map{|tf| experiments_by_tf(tf) }.select{|dataset| matrix.datasets.include?(dataset) }.uniq
  motifs = common_tfs.flat_map{|tf| motifs_by_tf(tf) }.select{|motif| matrix.motifs.include?(motif) }.uniq
  matrix.bound_at(datasets, motifs)
end

##############################

def dataset_groups_by_tf(datasets)
  datasets.group_by{|dataset| tf_by_experiment(dataset) }
end

def motif_groups_by_tf(motifs)
  result = Hash.new{|h,k| h[k] = []}
  motifs.each{|motif|
    tfs_by_motif(motif).each{|tf|
      result[tf] << motif
    }
  }
  result
end

def dataset_groups_by_family(datasets, aggregation_level)
  result = Hash.new{|h,k| h[k] = []}
  datasets.each{|dataset|
    tf = tf_by_experiment(dataset)
    families = tf_info_by_gene_name(tf)[aggregation_level]
    families.each{|family|
      result[family] << dataset
    }
  }
  result
end

def motif_groups_by_family(motifs, aggregation_level)
  result = Hash.new{|h,k| h[k] = [] }
  motifs.each{|motif|
    tfs_by_motif(motif).each{|tf|
      families = tf_info_by_gene_name(tf)[aggregation_level]
      families.each{|family|
        result[family] << motif
      }
    }
  }
  result.transform_values{|v| v.uniq }
end

##############################

def family_idx(family)
  family.match(/\{(\d+(\.\d+)*)\}/)[1].split('.').map(&:to_i)
end

def families_by_tf(tf, aggregation_level)
  tf_info_by_gene_name(tf)[aggregation_level].sort
end

def families_idx_by_tf(tf, aggregation_level)
  families_by_tf(tf, aggregation_level).map{|family| family_idx(family) }.sort
end

def families_by_dataset(dataset, aggregation_level)
  tf = tf_by_experiment(dataset)
  families_by_tf(tf, aggregation_level).sort
end

def families_idx_by_dataset(dataset, aggregation_level)
  families_by_dataset(dataset, aggregation_level).map{|family| family_idx(family) }.sort
end

def families_idx_by_motif(motif, aggregation_level)
  tfs_by_motif(motif).flat_map{|tf|
    tf_info_by_gene_name(tf)[aggregation_level]
  }.uniq.map{|fam| family_idx(fam) }.sort
end
