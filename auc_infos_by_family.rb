require 'json'
require_relative 'aucs'
require_relative 'tf_motifs_mapping'

class Array
  def mean; empty? ? nil : sum(0.0) / size; end
  def median
    return nil if empty?
    sorted = self.sort
    if size % 2 == 1
      sorted[size / 2]
    else
      (sorted[size / 2] + sorted[size / 2 - 1]) / 2.0
    end
  end
end

AUCS = Aucs.from_file('source_data/chipseq/motifs_vs_remap.tsv')

families = TF_INFO_BY_NAME.flat_map{|tf, tf_infos|
  tf_infos[:tf_family]
}.uniq

EXPERIMENTS_BY_FAMILY = families.map{|family|
  experiments = AUCS.experiments.select{|exp|
    tf = TF_BY_EXPERIMENT[exp]
    TF_INFO_BY_NAME[tf][:tf_family].include?(family)
  }
  [family, experiments]
}.to_h

MOTIFS_BY_FAMILY = families.map{|family|
  motifs = AUCS.motifs.select{|motif|
    tfs = TFS_BY_MOTIF[motif]
    tfs.flat_map{|tf|
      TF_INFO_BY_NAME[tf][:tf_family]
    }.include?(family)
  }
  [family, motifs]
}.to_h

#######################

representative_motifs = File.readlines('source_data/motifs/hocomoco_cluster_list_192.txt').map{|l|
  l.chomp.split("\t").first
}.reject{|motif|
  ['NFKB1_HUMAN.H11MO.0.A', 'EVX1_HUMAN.H11MO.0.D', 'EVX2_HUMAN.H11MO.0.A'].include?(motif)
}

REPRESENTATIVE_MOTIFS_BY_FAMILY = families.map{|family|
  motifs = representative_motifs.select{|motif|
    tfs = TFS_BY_MOTIF[motif] || []
    tfs.flat_map{|tf|
      TF_INFO_BY_NAME[tf][:tf_family]
    }.include?(family)
  }
  [family, motifs]
}.to_h

#######################

def aucs_over_family_datasets(motifs, family)
  experiments = EXPERIMENTS_BY_FAMILY[family] || []
  motifs.compact.flat_map{|motif|
    AUCS.motif_aucs_across_experiments(motif, experiments).values
  }
end

#######################

all_aucs_by_family = MOTIFS_BY_FAMILY.map{|family, motifs|
  [family, aucs_over_family_datasets(motifs, family)]
}.to_h

#######################

best_motif_by_family_in_family = families.map{|family|
  motif_aucs_pairs = MOTIFS_BY_FAMILY[family].map{|motif|
    [motif, aucs_over_family_datasets([motif], family)]
  }
  motif, auc_vals = motif_aucs_pairs.max_by{|motif, aucs| aucs.mean }
  [family, motif]
}.to_h

best_motif_aucs_by_family_in_family = families.map{|family|
  best_motif = best_motif_by_family_in_family[family]
  [family, aucs_over_family_datasets([best_motif], family)]
}.to_h

#######################

best_motif_by_family_overall = families.map{|family|
  motif_aucs_pairs = AUCS.motifs.map{|motif|
    [motif, aucs_over_family_datasets([motif], family)]
  }
  motif, auc_vals = motif_aucs_pairs.max_by{|motif, aucs| aucs.mean }
  [family, motif]
}.to_h

best_motif_aucs_by_family_overall = families.map{|family|
  best_motif = best_motif_by_family_overall[family]
  [family, aucs_over_family_datasets([best_motif], family)]
}.to_h

#######################

representative_motifs_aucs_by_family = families.map{|family|
  motifs_aucs = REPRESENTATIVE_MOTIFS_BY_FAMILY[family].map{|motif|
    [motif, aucs_over_family_datasets([motif], family)]
  }.to_h
  [family, motifs_aucs]
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
    representatives: representative_motifs_aucs_by_family[family],
  }
  [family, infos]
}.to_h

File.write('auc_infos_by_family.json', aucs_infos_by_family.to_json)


File.open('auc_infos_by_family.tsv', 'w') {|fw|
  fw.puts ['family', 'type', 'auc'].join("\t")
  aucs_infos_by_family.flat_map{|family, family_infos|
    rows = []
    rows += family_infos[:all].map{|auc| [family, 'all', auc] }
    rows += family_infos[:best_in_family].map{|auc| [family, 'best_in_family', auc] }
    rows += family_infos[:best_overall].map{|auc| [family, 'best_overall', auc] }
    rows += family_infos[:representatives].flat_map{|motif, motif_aucs|
      motif_aucs.map{|auc|
        [family, "representative:#{motif}", auc]
      }
    }
    rows += family_infos[:representatives].flat_map{|motif, motif_aucs|
      motif_aucs.map{|auc|
        [family, "all_representatives", auc]
      }
    }
    rows
  }.each{|row|
    fw.puts row.join("\t")
  }
}

#######################

auc_delta_infos = families.map{|family|
  best_motif_overall = best_motif_by_family_overall[family]
  best_motif_in_family = best_motif_by_family_in_family[family]
  representative_motifs = REPRESENTATIVE_MOTIFS_BY_FAMILY[family]

  if best_motif_in_family
    best_in_family_delta_aucs = EXPERIMENTS_BY_FAMILY[family].map{|experiment|
      AUCS.auc(best_motif_overall, experiment) - AUCS.auc(best_motif_in_family, experiment)
    }
  else
    best_in_family_delta_aucs = []
  end

  representatives_infos = representative_motifs.map{|representative_motif|
    representative_delta_aucs = EXPERIMENTS_BY_FAMILY[family].map{|experiment|
      AUCS.auc(best_motif_overall, experiment) - AUCS.auc(representative_motif, experiment)
    }
    [representative_motif, representative_delta_aucs]
  }.to_h


  [family, {best_in_family_delta_aucs: best_in_family_delta_aucs, representatives_infos: representatives_infos}]
}.to_h

#######################

File.write('auc_delta_infos.json', auc_delta_infos.to_json)

File.open('auc_delta_infos.tsv', 'w') {|fw|
  fw.puts ['family', 'delta_type', 'auc'].join("\t")
  auc_delta_infos.flat_map{|family, family_infos|
    rows = []
    rows += family_infos[:best_in_family_delta_aucs].map{|auc|
      [family, 'best_in_family_delta', auc]
    }
    rows += family_infos[:representatives_infos].sort_by{|motif,repr_aucs|
      repr_aucs.mean
    }.each_with_index.flat_map{|(motif, repr_aucs), idx|
      repr_aucs.map{|auc|
        [family, "representative_delta:#{idx}", auc]
      }
    }
    rows
  }.each{|row|
    fw.puts row.join("\t")
  }
}
