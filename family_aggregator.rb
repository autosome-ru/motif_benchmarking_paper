RETAIN = [
  'unknown',
  # 'ambiguous',

  # TFClass classes
  'Basic leucine zipper factors (bZIP){1.1}',
  'C2H2 zinc finger factors{2.3}',
  'Tryptophan cluster factors{3.5}',
  'Basic helix-loop-helix factors (bHLH){1.2}',
  'Homeo domain factors{3.1}',
  'Nuclear receptors with C4 zinc fingers{2.1}',
  'Fork head / winged helix factors{3.3}',
  'Other C4 zinc finger-type factors{2.2}',
  'High-mobility group (HMG) domain factors{4.1}',
  'Uncharacterized{0.0}',
  'Heteromeric CCAAT-binding factors{4.2}',
  'Rel homology region (RHR) factors{6.1}',
  'Runt domain factors{6.4}',
  'MADS box factors{5.1}',
  'STAT domain factors{6.2}',
  'SMAD/NF-1 DNA-binding domain factors{7.1}',

  # TFClass families
  'Ets-related factors{3.5.2}',
  'More than 3 adjacent zinc finger factors{2.3.3}',
  'NRF{0.0.6}',
  'Three-zinc finger Krüppel-related factors{2.3.1}',
  'GATA-type zinc fingers{2.2.1}',
  'Tal-related factors{1.2.3}',
  'bHLH-ZIP factors{1.2.6}',
  'Forkhead box (FOX) factors{3.3.1}',
  'More than 3 adjacent zinc finger factors{2.3.3}',
]

def retain_or_other(family)
  RETAIN.include?(family) ? family : 'other'
end

drop_unknown_experiment = ARGV.delete('--drop-unknown-experiment')
drop_unknown_motif = ARGV.delete('--drop-unknown-motif')
filename = ARGV[0]

header = ['experiment_TF_family', 'best_motifs_family']
puts header.join("\t")
File.readlines(filename).map{|l|
  experiment_TF, experiment_TF_family, best_motif, best_auc, tfs_of_best_motif, best_motifs_family = l.chomp.split("\t")
  [experiment_TF_family, best_motifs_family]
}.reject{|experiment_TF_family, best_motifs_family|
  drop_unknown_experiment && experiment_TF_family == 'unknown'
}.reject{|experiment_TF_family, best_motifs_family|
  drop_unknown_motif && best_motifs_family == 'unknown'
}.map{|experiment_TF_family, best_motifs_family|
  [experiment_TF_family, best_motifs_family].map{|fam|
    retain_or_other(fam)
  }
}.each{|row|
  puts row.join("\t")
}
