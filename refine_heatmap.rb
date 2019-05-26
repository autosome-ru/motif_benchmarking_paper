require 'json'
require_relative 'annotation'
require_relative 'aucs_matrix'

use_classnames = ARGV.delete('--class-names')
heatmap_fn = ARGV[0]

aucs_matrix_chipseq = AucsMatrix.from_file('source_data/chipseq/motifs_vs_remap.tsv')
chipseq_annotation = Annotation.new(aucs_matrix_chipseq.experiments, aucs_matrix_chipseq.motifs)

aucs_matrix_selex = AucsMatrix.from_file('source_data/selex/motifs_vs_selex10.tsv')
selex_annotation = Annotation.new(aucs_matrix_selex.experiments, aucs_matrix_selex.motifs)

lines = File.readlines(heatmap_fn).map(&:chomp)
families = lines.first.split("\t").drop(1)

matrix = lines.drop(1).map{|l|
  l.split("\t", 1 + families.size).drop(1).map{|x|
    x.empty? ? nil : Float(x)
  }
}

relevant_families = chipseq_annotation.tfclass_names(:tf_family).select{|family|
  chipseq_annotation.motifs_by_tfclass_name(:tf_family, family).size >= 2
}.select{|family|
  selex_annotation.experiments_by_tfclass_name(:tf_family, family).size >= 2
}.select{|family|
  chipseq_annotation.experiments_by_tfclass_name(:tf_family, family).size >= 2
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
