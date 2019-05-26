require_relative 'aucs_matrix'
require_relative 'auc_collection'

class Array
  def mean; empty? ? nil : sum(0.0) / size; end
end

drop_empty = ARGV.delete('--drop-empty')
aucs_matrix = ARGV[0]
tfclass_level = ARGV[1].to_sym # :tf_class

aucs_matrix = AucsMatrix.from_file(aucs_matrix)
annotation = Annotation.new(aucs_matrix.experiments, aucs_matrix.motifs)
aucs = AucCollection.new(aucs_matrix, annotation)


families = annotation.tfclass_names(tfclass_level)

tf_classes = families
# .map{|family|
#   annotation.tfclass_parent(family)
# }.chunk(&:itself).map{|group_name, vs|
#   [group_name] + ['.'] * (vs.size - 1)
# }.flatten

header = ['-', *tf_classes]
heatmap = [header] + families.zip(tf_classes).map{|motif_family, motif_tfclass|
  motifs = annotation.motifs_by_tfclass_name(tfclass_level, motif_family)
  row = families.map{|experiment_family|
    experiments = annotation.experiments_by_tfclass_name(tfclass_level, experiment_family)
    aucs.aucs_subset(motifs, experiments).mean
  }
  [motif_tfclass, *row]
}

# Drop empty rows/columns in a matrix
if drop_empty
  heatmap = heatmap.select{|motif_family, *row|
    row.any?
  }.transpose.select{|experiment_family, *column|
    column.any?
  }.transpose
end

heatmap.each{|row|
  puts row.join("\t")
}
