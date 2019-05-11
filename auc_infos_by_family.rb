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

aucs = Aucs.from_file('source_data/chipseq/motifs_vs_remap.tsv')

experiments_by_family = aucs.experiments.select{|exp|
  tf = TF_BY_EXPERIMENT[exp]
  TF_INFO_BY_NAME[tf][:tf_family].size == 1
}.group_by{|exp|
  tf = TF_BY_EXPERIMENT[exp]
  TF_INFO_BY_NAME[tf][:tf_family].first
}

motifs_by_family = aucs.motifs.reject{|motif|
  TFS_BY_MOTIF[motif].size != 1
}.group_by{|motif|
  tf = TFS_BY_MOTIF[motif].first
  TF_INFO_BY_NAME[tf][:tf_family].first
}.to_h

aucs_by_family = motifs_by_family.map{|family, motifs|
  experiments = experiments_by_family[family] || []
  auc_vals = motifs.flat_map{|motif|
    experiments.map{|experiment|
      aucs.auc(motif, experiment)
    }
  }
  [family, auc_vals]
}.to_h

best_motif_by_family = motifs_by_family.map{|family, motifs|
  experiments = experiments_by_family[family] || []
  motif, auc_vals = motifs.map{|motif|
    vals = experiments.map{|experiment|
      aucs.auc(motif, experiment)
    }
    [motif, vals]
  }.max_by{|motif, vals|
    vals.mean
  }
  [family, motif]
}.to_h

best_motif_aucs_by_family = best_motif_by_family.map{|family, motif|
  experiments = experiments_by_family[family] || []
  motif_aucs = experiments.map{|experiment|
    aucs.auc(motif, experiment)
  }
  [family, motif_aucs]
}.to_h

representative_motifs = File.readlines('source_data/motifs/hocomoco_cluster_list_192.txt').map{|l|
  l.chomp.split("\t").first
}.reject{|motif|
  ['NFKB1_HUMAN.H11MO.0.A', 'EVX1_HUMAN.H11MO.0.D', 'EVX2_HUMAN.H11MO.0.A'].include?(motif)
}

representative_motifs_by_family = representative_motifs.select{|motif|
  tfs = TFS_BY_MOTIF[motif] || []
  fams = tfs.flat_map{|tf|
    TF_INFO_BY_NAME[tf][:tf_family]
  }.uniq
  fams.size == 1
}.group_by{|motif|
  tfs = TFS_BY_MOTIF[motif] || []
  fams = tfs.flat_map{|tf|
    TF_INFO_BY_NAME[tf][:tf_family]
  }.uniq
  fams.first
}

representative_motifs_aucs_by_family = representative_motifs_by_family.map{|family, motifs|
  experiments = experiments_by_family[family] || []
  motifs_aucs = motifs.map{|motif|
    motif_aucs = experiments.map{|experiment|
      aucs.auc(motif, experiment)
    }
    [motif, motif_aucs]
  }.to_h
  [family, motifs_aucs]
}.to_h

aucs_infos_by_family = [
  'Ets-related factors{3.5.2}',
  'Three-zinc finger Krüppel-related factors{2.3.1}',
  'Factors with multiple dispersed zinc fingers{2.3.4}',
  'Forkhead box (FOX) factors{3.3.1}',
].map{|family|
  infos = {
    all: aucs_by_family[family],
    best: best_motif_aucs_by_family[family],
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
    rows += family_infos[:best].map{|auc| [family, 'best', auc] }
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
