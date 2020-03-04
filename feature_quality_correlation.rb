def normalize_motif_name(motif_name)
  if motif_name.match(/\.H11MO.\d/) # hocomoco
    motif_name
  elsif motif_name.start_with?('MA') # jaspar
    motif_name.split('_')[0]
  else # cisbp
    motif_name.split('_')[0,2].join('_')
  end
end

# filename = ARGV[0] # 'remap_motif_tf.tsv'

motif_infos = File.readlines('source_data/annotation/motif_features_uniq.tsv').drop(1).map{|l|
  collection, motif, source_type, motif_length, mean_ic, gc_content = l.chomp.split("\t")
  [normalize_motif_name(motif), {collection: collection, source_type: source_type, motif_length: motif_length.to_i, mean_ic: mean_ic.to_f, gc_content: gc_content.to_f}]
}.to_h

header = ['auc', 'motif_length', 'gc_content', 'mean_ic']
puts header.join("\t")

input_header = readline 
ARGF.each_line.map{|l|
  motif, auc, experiment_tf, motif_tfs = l.chomp.split("\t")
  [normalize_motif_name(motif), auc.to_f]
}.select{|motif, auc|
  motif_infos.has_key?(motif)
}.each{|motif, auc|
  infos = [auc, motif_infos[motif].values_at(:motif_length, :gc_content, :mean_ic)]
  puts infos.join("\t")
}
