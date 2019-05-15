require 'json'

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

TFS_BY_MOTIF = motif_tfs.to_h
MOTIFS_BY_TF = tf_motifs.to_h

#####################

tf_infos = File.readlines('all_tf_infos.json').map{|l|
  JSON.parse(l, symbolize_names: true)
}

TF_INFO_BY_NAME = tf_infos.group_by{|info|
  info[:tf_gene_name]
}.map{|k,vs|
  [k, vs.first]
}.to_h

#####################

experiment_tf_pairs = [
  'source_data/chipseq/remap_genes2exp.txt',
  'source_data/selex/jolma13_genes2exp.txt',
].flat_map{|experiment_mapping_fn|
  File.readlines(experiment_mapping_fn).map{|l|
    tf, experiment = l.chomp.split("\t").first(2)
    [experiment, tf]
  }
}

TF_BY_EXPERIMENT = experiment_tf_pairs.to_h
EXPERIMENTS_BY_TF = experiment_tf_pairs.group_by{|experiment, tf|
  tf
}.map{|tf, exp_tf_pairs|
  [tf, exp_tf_pairs.map{|exp,tf| exp }]
}.to_h
