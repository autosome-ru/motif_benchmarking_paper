require 'json'
require_relative 'annotation'
require_relative 'aucs_matrix'
require_relative 'auc_collection'

class Array
  def mean; empty? ? nil : sum(0.0) / size; end
end

def certain_or_ambiguous(values)
  if values.size == 1
    values.first
  elsif values.size > 1
    'ambiguous'
  elsif values.size == 0
    'unknown'
  end
end



aucs_matrix = AucsMatrix.from_file('source_data/chipseq/motifs_vs_remap.tsv')
annotation = Annotation.new(aucs_matrix.experiments, aucs_matrix.motifs)
aucs = AucCollection.new(aucs_matrix, annotation)

performance_infos = annotation.motif_tfs.flat_map{|motif_tf|
  motifs = annotation.motifs_by_tf(motif_tf)
  annotation.experiment_tfs.map{|experiment_tf|
    experiments = annotation.experiments_by_tf(experiment_tf)
    
    best_motif, best_mean_auc = motifs.map{|motif|
      [motif, aucs.aucs_subset([motif], experiments).mean]
    }.max_by{|motif, mean_auc|
      mean_auc
    }

    {motif_tf: motif_tf, experiment_tf: experiment_tf, motif: best_motif, auc: best_mean_auc}
  }
}

performance_infos = performance_infos.map{|info|
  motif_tf_info = annotation.tf_info_by_gene_name(info[:motif_tf])
  motif_tf_family = certain_or_ambiguous(motif_tf_info[:tf_family])
  motif_tf_class = certain_or_ambiguous(motif_tf_info[:tf_class])

  experiment_tf_info = annotation.tf_info_by_gene_name(info[:experiment_tf])
  experiment_tf_family = certain_or_ambiguous(experiment_tf_info[:tf_family])
  experiment_tf_class = certain_or_ambiguous(experiment_tf_info[:tf_class])

  info.merge({
    motif_family: motif_tf_family,
    experiment_family: experiment_tf_family,
    motif_class: motif_tf_class,
    experiment_class: experiment_tf_class,
  })
}

column_order = [:motif, :motif_tf, :motif_family, :motif_class, :experiment_tf, :experiment_family, :experiment_class, :auc]

File.open('meanauc_0.8.tsv', 'w'){|fw|
  fw.puts column_order.join("\t")
  performance_infos.select{|info|
    info[:auc] >= 0.8
  }.each{|infos|
    fw.puts infos.values_at(*column_order).join("\t")
  }
}


# annotation.motifs.flat_map{|motif|
#   annotation.experiment_tfs.map{|experiment_tf|
#     experiments = annotation.experiments_by_tf(experiment_tf)
#     mean_auc = aucs.aucs_subset([motif], experiments).mean
#     [experiment_tf, motif, mean_auc]
#   }
# }

# .each{|row|
#   puts row.join("\t")
# }
