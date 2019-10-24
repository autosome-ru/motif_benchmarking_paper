require 'json'
require_relative 'aucs_matrix'
require_relative 'auc_collection'

class Array
  def mean; empty? ? nil : sum(0.0) / size; end
end

aucs_matrix = AucsMatrix.from_file('source_data/chipseq/motifs_vs_remap.tsv')
annotation = Annotation.new(aucs_matrix.experiments, aucs_matrix.motifs)
aucs = AucCollection.new(aucs_matrix, annotation)

#######################
families = annotation.tfclass_names(:tf_family)

representative_motifs = File.readlines('source_data/motifs/representatives.txt').map(&:chomp)

REPRESENTATIVE_MOTIFS_BY_FAMILY = families.map{|family|
  motifs = representative_motifs.select{|motif|
    tfs = annotation.tfs_by_motif(motif)
    tfs.flat_map{|tf|
      annotation.tf_info_by_gene_name(tf)[:tf_family]
    }.include?(family)
  }
  [family, motifs]
}.to_h

def representative_motifs_by_family(aucs, family)
  @repr_mots_by_fam ||= REPRESENTATIVE_MOTIFS_BY_FAMILY.map{|fam, motifs|
    [fam, motifs.select{|motif| aucs.motifs.include?(motif) }]
  }.to_h
  @repr_mots_by_fam[family]
end

#######################


family_triples = families.map{|family|
  motifs = annotation.motifs_by_tfclass_name(:tf_family, family)
  experiments = annotation.experiments_by_tfclass_name(:tf_family, family)
  [family, motifs, experiments]
}

all_aucs_by_family = family_triples.map{|family, motifs, experiments|
  [family, aucs.aucs_subset(motifs, experiments)]
}.to_h

#######################

best_motif_by_family_in_family = family_triples.map{|family, motifs, experiments|
  motif_aucs_pairs = motifs.map{|motif|
    [motif, aucs.aucs_subset([motif], experiments)]
  }
  motif, auc_vals = motif_aucs_pairs.max_by{|motif, aucs| aucs.mean }
  [family, motif]
}.to_h

best_motif_aucs_by_family_in_family = family_triples.map{|family, motifs, experiments|
    best_motif = best_motif_by_family_in_family[family]
    [family, aucs.aucs_subset([best_motif].compact, experiments)]
}.to_h

#######################

best_motif_by_family_overall = family_triples.map{|family, _, experiments|
  motif_aucs_pairs = annotation.motifs.map{|motif|
    [motif, aucs.aucs_subset([motif], experiments)]
  }
  motif, auc_vals = motif_aucs_pairs.max_by{|motif, aucs| aucs.mean }
  [family, motif]
}.to_h

best_motif_aucs_by_family_overall = family_triples.map{|family, _, experiments|
  best_motif = best_motif_by_family_overall[family]
  [family, aucs.aucs_subset([best_motif].compact, experiments)]
}.to_h

#######################

best_representative_by_family = family_triples.map{|family, _, experiments|
  motif_aucs_pairs = representative_motifs.map{|motif|
    [motif, aucs.aucs_subset([motif], experiments)]
  }
  motif, auc_vals = motif_aucs_pairs.max_by{|motif, aucs| aucs.mean }
  [family, motif]
}.to_h

best_representative_aucs_by_family = family_triples.map{|family, _, experiments|
  best_motif = best_representative_by_family[family]
  [family, aucs.aucs_subset([best_motif], experiments)]
}.to_h

#######################

family_representative_motifs_aucs_by_family = family_triples.map{|family, _, experiments|
  motifs_aucs = representative_motifs_by_family(aucs, family).map{|motif|
    [motif, aucs.aucs_subset([motif], experiments)]
  }.to_h
  [family, motifs_aucs]
}.to_h

representative_motifs_aucs_by_family = family_triples.map{|family, _, experiments|
  motifs_aucs = representative_motifs.map{|motif|
    [motif, aucs.aucs_subset([motif], experiments)]
  }.to_h
  [family, motifs_aucs]
}.to_h


#######################

proper_tf_aucs_by_family = families.map{|family|
  tfs = annotation.tfs_by_tfclass_name(:tf_family, family)
  tf_aucs = tfs.map{|tf|
    [tf, aucs.aucs_subset(annotation.motifs_by_tf(tf), annotation.experiments_by_tf(tf))]
  }.to_h
  [family, tf_aucs]
}.to_h

#######################

selected_families = [
  'Ets-related factors{3.5.2}',
  'Three-zinc finger Krüppel-related factors{2.3.1}',
  'Factors with multiple dispersed zinc fingers{2.3.4}',
  'Forkhead box (FOX) factors{3.3.1}',
]

aucs_infos_by_family = selected_families.map{|family|
  infos = {
    all: all_aucs_by_family[family],
    best_in_family: best_motif_aucs_by_family_in_family[family],
    best_overall: best_motif_aucs_by_family_overall[family],
    best_representative: best_representative_aucs_by_family[family],
    family_representatives: family_representative_motifs_aucs_by_family[family],
    proper_tf: proper_tf_aucs_by_family[family],
  }
  [family, infos]
}.to_h

# File.write('auc_infos_by_family.json', aucs_infos_by_family.to_json)


File.open('auc_infos_by_family.tsv', 'w') {|fw|
  fw.puts ['family', 'type', 'auc'].join("\t")
  aucs_infos_by_family.flat_map{|family, family_infos|
    rows = []
    rows += family_infos[:all].map{|auc| [family, 'all', auc] }
    rows += family_infos[:best_in_family].map{|auc| [family, 'best_in_family', auc] }
    rows += family_infos[:best_overall].map{|auc| [family, 'best_overall', auc] }
    rows += family_infos[:best_representative].map{|auc| [family, 'best_representative', auc] }
    rows += family_infos[:proper_tf].flat_map{|tf, aucs|
      aucs.flat_map{|auc|
        [
          [family, "proper_tf:#{tf}", auc],
          [family, "proper_tfs", auc],
        ]
      }
    }

    rows += family_infos[:family_representatives].sort_by{|motif, motif_aucs|
      motif_aucs.mean
    }.reverse.each_with_index.flat_map{|(motif, motif_aucs), idx|
      motif_aucs.map{|auc|
        [family, "representative:#{idx}:#{motif}", auc]
      }
    }
    rows += family_infos[:family_representatives].flat_map{|motif, motif_aucs|
      motif_aucs.map{|auc|
        [family, "all_representatives", auc]
      }
    }
    rows
  }.each{|row|
    fw.puts row.join("\t")
  }
}

# #######################

# auc_delta_infos = family_triples.map{|family, motifs, experiments|
#   best_motif_overall = best_motif_by_family_overall[family]
#   best_motif_in_family = best_motif_by_family_in_family[family]
#   family_representative_motifs = representative_motifs_by_family(aucs, family)

#   if best_motif_in_family
#     best_in_family_delta_aucs = experiments.map{|experiment|
#       aucs.auc(best_motif_overall, experiment) - aucs.auc(best_motif_in_family, experiment)
#     }
#   else
#     best_in_family_delta_aucs = []
#   end

#   representatives_infos = family_representative_motifs.map{|representative_motif|
#     representative_delta_aucs = experiments.map{|experiment|
#       aucs.auc(best_motif_overall, experiment) - aucs.auc(representative_motif, experiment)
#     }
#     [representative_motif, representative_delta_aucs]
#   }.to_h


#   [family, {best_in_family_delta_aucs: best_in_family_delta_aucs, representatives_infos: representatives_infos}]
# }.to_h

# #######################

# File.write('auc_delta_infos.json', auc_delta_infos.to_json)

# File.open('auc_delta_infos.tsv', 'w') {|fw|
#   fw.puts ['family', 'delta_type', 'auc'].join("\t")
#   auc_delta_infos.flat_map{|family, family_infos|
#     rows = []
#     rows += family_infos[:best_in_family_delta_aucs].map{|auc|
#       [family, 'best_in_family_delta', auc]
#     }
#     rows += family_infos[:representatives_infos].sort_by{|motif,repr_aucs|
#       repr_aucs.mean
#     }.each_with_index.flat_map{|(motif, repr_aucs), idx|
#       repr_aucs.map{|auc|
#         [family, "representative_delta:#{idx}", auc]
#       }
#     }
#     rows
#   }.each{|row|
#     fw.puts row.join("\t")
#   }
# }
