require 'json'
require_relative 'aucs_matrix'
require_relative 'auc_collection'

class Array
  def mean; empty? ? nil : sum(0.0) / size; end
end

aucs_matrix_chipseq = AucsMatrix.from_file('source_data/chipseq/motifs_vs_remap.tsv');nil
chipseq_annotation = Annotation.new(aucs_matrix_chipseq.experiments, aucs_matrix_chipseq.motifs);nil
aucs_chipseq = AucCollection.new(aucs_matrix_chipseq, chipseq_annotation);nil

aucs_matrix_selex = AucsMatrix.from_file('source_data/selex/motifs_vs_selex10.tsv');nil
selex_annotation = Annotation.new(aucs_matrix_selex.experiments, aucs_matrix_selex.motifs);nil
aucs_selex = AucCollection.new(aucs_matrix_selex, selex_annotation);nil

tfs_w_chipseq_and_selex_motifs = chipseq_annotation.motif_tfs.select{|tf|
  types = chipseq_annotation.motifs_by_tf(tf).map{|motif| chipseq_annotation.motif_source_type(motif) }.uniq
  types.include?('ChIP-seq') && (types.include?('HT-SELEX') || types.include?('SMiLE-seq') || types.include?('SELEX') )
}

tfs_to_analyze = tfs_w_chipseq_and_selex_motifs.reject{|tf|
  chipseq_annotation.experiments_by_tf(tf).empty? || selex_annotation.experiments_by_tf(tf).empty?
}

columns = [:tf, 
  :num_chipseq_motifs, :num_selex_motifs, :num_chipseq_experiments, :num_selex_experiments,
  :auc_chipseq_on_chipseq, :auc_chipseq_on_selex, :auc_selex_on_chipseq, :auc_selex_on_selex,
]

puts columns.join("\t")
tfs_to_analyze.map{|tf|
  motifs = chipseq_annotation.motifs_by_tf(tf)
  chipseq_motifs = motifs.select{|motif|
    chipseq_annotation.motif_source_type(motif) == 'ChIP-seq'
  }
  selex_motifs = motifs.select{|motif|
    ['HT-SELEX', 'SELEX', 'SMiLE-seq'].include?(chipseq_annotation.motif_source_type(motif))
  }
  chipseq_experiments = chipseq_annotation.experiments_by_tf(tf)
  selex_experiments = selex_annotation.experiments_by_tf(tf)
  # p({tf: tf, chipseq_motifs: chipseq_motifs, selex_motifs: selex_motifs})
  { tf: tf,
    num_chipseq_motifs: chipseq_motifs.size,
    num_selex_motifs: selex_motifs.size,
    num_chipseq_experiments: chipseq_experiments.size,
    num_selex_experiments: selex_experiments.size,
    auc_chipseq_on_chipseq: aucs_chipseq.aucs_subset(chipseq_motifs, chipseq_experiments).mean,
    auc_chipseq_on_selex: aucs_selex.aucs_subset(chipseq_motifs, selex_experiments).mean,
    auc_selex_on_chipseq: aucs_chipseq.aucs_subset(selex_motifs, chipseq_experiments).mean,
    auc_selex_on_selex: aucs_selex.aucs_subset(selex_motifs, selex_experiments).mean,
  }
}.each{|infos|
  puts infos.values_at(*columns).join("\t")
}
