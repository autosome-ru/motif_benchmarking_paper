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

def motif_tf(motif)
  TF_BY_MOTIF[motif]
end
