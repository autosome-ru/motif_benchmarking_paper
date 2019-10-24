require 'json'
require_relative 'annotation'
require_relative 'aucs_matrix'

use_classnames = ARGV.delete('--class-names')
heatmap_fn = ARGV[0]
family_field = ARGV[1].to_sym # :tf_family / :cisbp_families

aucs_matrix_chipseq = AucsMatrix.from_file('source_data/final/remap_all_roc.txt')
chipseq_annotation = Annotation.new(aucs_matrix_chipseq.experiments, aucs_matrix_chipseq.motifs)

aucs_matrix_selex = AucsMatrix.from_file('source_data/final/jy50_all_roc.txt')
selex_annotation = Annotation.new(aucs_matrix_selex.experiments, aucs_matrix_selex.motifs)

lines = File.readlines(heatmap_fn).map(&:chomp)
families = lines.first.split("\t").drop(1)

matrix = lines.drop(1).map{|l|
  l.split("\t", 1 + families.size).drop(1).map{|x|
    x.empty? ? nil : Float(x)
  }
}

relevant_families = chipseq_annotation.tfclass_names(family_field).select{|family|
  chipseq_annotation.motifs_by_tfclass_name(family_field, family).size >= 2
}.select{|family|
  selex_annotation.experiments_by_tfclass_name(family_field, family).size >= 2
}.select{|family|
  chipseq_annotation.experiments_by_tfclass_name(family_field, family).size >= 2
}

relevant_fam_indices = families.each_with_index.select{|family,idx| relevant_families.include?(family) }.map(&:last)
families = families.values_at(*relevant_fam_indices)

matrix = matrix.values_at(*relevant_fam_indices).map{|row|
  row.values_at(*relevant_fam_indices)
}

if use_classnames
  family_labels = families.map{|family|
    chipseq_annotation.tfclass_parent(family)
  }.chunk(&:itself).map{|group_name, vs|
    [group_name] + ['--skip-me--'] * (vs.size - 1)
  }.flatten
else
  family_labels = families
end

puts ['-', *family_labels].join("\t")
matrix.zip(family_labels).map{|row, family|
  puts [family, *row].join("\t")
}
