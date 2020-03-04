require 'set'
require_relative 'aucs_matrix'
require_relative 'annotation'

def normalize_motif_name(motif_name)
  if motif_name.match(/\.H11MO.\d/) # hocomoco
    motif_name
  elsif motif_name.start_with?('MA') # jaspar
    motif_name.split('_')[0]
  else # cisbp
    motif_name.split('_')[0,2].join('_')
  end
end

aucs_matrix_fn = 'source_data/final/remap_all_roc.txt'

aucs_matrix = AucsMatrix.from_file(aucs_matrix_fn); nil

motifs = aucs_matrix.motifs.map{|motif| normalize_motif_name(motif) }.to_set
File.open('source_data/annotation/motif_annotation_final_features_cleanup.tsv', 'w'){|fw|
  File.open('source_data/annotation/motif_annotation_final_features.tsv'){|f|
    header = f.readline
    fw.puts(header)
    f.each_line.map{|l|
      motif = normalize_motif_name(l.chomp.split("\t")[1])
      fw.puts(l) if motifs.include?(motif)
    }
  }
};nil
