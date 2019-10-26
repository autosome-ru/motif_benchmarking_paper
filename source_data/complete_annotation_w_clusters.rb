motif_centers = File.readlines('annotation/hocomoco_cluster_list_192.txt').flat_map{|l|
  center, cluster_els = l.chomp.split("\t")
  cluster_els.split(';').map{|el|
    [el, center]
  }
}.to_h


additional_motif_centers = File.readlines('annotation/missing_cluster_representatives.tsv').map{|l|
  motif, cluster_representative, sim = l.chomp.split("\t").map(&:strip)
  [motif, cluster_representative, sim.to_f]
}.select{|motif, cluster_representative, sim|
  sim >= 0.05
}.map{|motif, cluster_representative, sim|
  center = motif_centers[cluster_representative]
  [motif, center]
}.to_h

File.readlines('annotation/hoco_ann.tsv').each{|l|
  motif, tf, family, subfamily, center, is_central = l.chomp.split("\t").map(&:strip)
  center = additional_motif_centers.fetch(motif, '-')  if center == '-'
  is_central = false
  if center != '-'
    cmd = "java -cp ape-3.0.2.jar ru.autosome.macroape.EvalSimilarity all_hoco_motifs/#{motif}.pwm all_hoco_motifs/#{center}.pwm"
    similarity = `#{cmd}`.lines.detect{|l| l.start_with?("S\t") }.split("\t").last.to_f
    is_central = (similarity == 1.0) # sometimes different motifs of a cluster (human and mouse) are the same motif
  end
  
  infos = [motif, tf, family, subfamily, center, is_central ? 'center' : 'periphery']
  puts infos.join("\t")
}
