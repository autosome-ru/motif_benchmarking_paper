motif_centers = File.readlines('annotation/hocomoco_cluster_list_192.txt').flat_map{|l|
  center, cluster_els = l.chomp.split("\t")
  cluster_els.split(';').map{|el|
    [el, center]
  }
}.to_h

ARGF.each_line{|l|
  l.chomp!
  motif = l.split("\t").first
  center = motif_centers.fetch(motif, '-')
  centrality = (motif == center) ? 'center' : 'periphery'
  puts "#{l}\t#{center}\t#{centrality}"
}
