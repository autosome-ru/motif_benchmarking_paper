require 'json'

tfclass_level_name = ARGV[0].to_sym

motif_tf_pairs ||= [
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

motif_tfs = motif_tf_pairs.group_by{|tf, motif|
  motif
}.map{|motif, pairs|
  [motif, pairs.map{|tf, motif| tf }]
}

tf_infos = File.readlines('all_tf_infos.json').map{|l|
  JSON.parse(l, symbolize_names: true)
}

tf_info_by_name = tf_infos.group_by{|info|
  info[:tf_gene_name]
}.map{|k,vs|
  [k, vs.first]
}.to_h

motif_families = motif_tfs.map{|motif, tfs|
  fams = tfs.map{|tf|
    tf_info_by_name[tf]
  }.flat_map{|info|
    info[tfclass_level_name]
  }.uniq
  [motif, fams]
}

motif_families.each{|motif, fams|
  if fams.size == 1
    fam = fams.first
  elsif fams.size == 0
    fam = 'unknown'
  elsif fams.size > 1
    fam = 'ambiguous'
  end
  puts [motif, fam].join("\t")
}
