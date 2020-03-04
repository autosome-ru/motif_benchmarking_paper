Signal.trap('PIPE', 'EXIT')
require 'json'
require_relative 'aucs_matrix'
require_relative 'annotation'

class Array
  def mean; empty? ? nil : sum(0.0) / size; end
  def median; empty? ? nil : sort[size/2]; end
end

# aucs_matrix_fn = 'source_data/final/remap_all_roc.txt'
aucs_matrix_fn = 'source_data/final/jy10_all_roc.txt'
tfclass_level_name = :tf_subfamily

aucs_matrix = AucsMatrix.from_file(aucs_matrix_fn)
annotation = Annotation.new(aucs_matrix.experiments, aucs_matrix.motifs)

experiment_tfs = annotation.experiment_tfs

infos = annotation.motifs.flat_map{|motif|
  tfs = annotation.tfs_by_motif(motif)
  tfs.flat_map{|tf|
    tfclass_names = annotation.tf_info_by_gene_name(tf)[tfclass_level_name]
    tf_experiments = annotation.experiments_by_tf(tf)
    tf_aucs = tf_experiments.map{|experiment|
      aucs_matrix.auc(motif, experiment)
    }
    tfclass_names.map{|tfclass_name|
      family_experiments = annotation.experiments_by_tfclass_name(tfclass_level_name, tfclass_name)
      family_aucs = family_experiments.map{|experiment|
        aucs_matrix.auc(motif, experiment)
      }
      [motif, tf,  tf_aucs.mean, tf_aucs.median, tfclass_name, family_aucs.mean, family_aucs.median]
    }
  }
}
# .select{|motif, tf, tf_mean, tf_median, family, family_mean, family_median|
#   family_mean
# }

# infos.group_by{|motif, tf, tf_mean, tf_median, family, family_mean, family_median|
#   family
# }.map{|family, rows|
#   rows.sort_by{|motif, tf, tf_mean, tf_median, family, family_mean, family_median| -family_mean }.first
# }.each{|row|
#   puts row.join("\t")
# }
header = ['motif', 'tf', 'tf_mean', 'tf_median', 'family', 'family_mean', 'family_median']
puts header.join("\t")
infos.each{|row|
  puts row.join("\t")
}