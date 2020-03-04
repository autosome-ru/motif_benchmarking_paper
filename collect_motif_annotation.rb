require 'json'
require 'httparty'

def request_uniprot_data_by_uniprot_ac(uniprot_ac)
  if !$cache
    if File.exist?('uniprot_cache.json')
      $cache = JSON.parse(File.read('uniprot_cache.json')).transform_values{|v| v.transform_keys(&:to_sym) }.to_h
    else
      $cache = {}
    end
  end
  return $cache[uniprot_ac]  if $cache[uniprot_ac]
  $cache[uniprot_ac] ||= begin
    base_url = "https://www.uniprot.org/uniprot/"
    params = {
      # organism:9606 - Homo sapiens (Human); Text query for organism obtains genes of other species 
      # like "Human adenovirus" organism on par with Human itself
      query: "accession:#{uniprot_ac}",
      sort: 'score',
      columns: 'id,entry name,genes',
      # limit: 1,
      format: 'tab'
    }
    HTTParty.get(base_url, query: params).body.lines.drop(1).map{|l|
      uniprot_ac, uniprot_id, gene_names = l.chomp.split("\t", 3).map(&:strip)
      gene_names = gene_names.split(/\s+/).map(&:strip)
      {uniprot_ac: uniprot_ac, uniprot_id: uniprot_id, gene_names: gene_names}
    }.first
  end
  File.write('uniprot_cache.json', $cache.to_json)
  $cache[uniprot_ac]
end

tf_to_cisbp = File.readlines('source_data/annotation/cisbp_refined.tsv').flat_map{|l|
  motif, motif_name, gene, cisbp_families, cisbp_dbds, cisbp_uniprot_acs, family_1, family_2, family_3, family_4, source_type, sr_model, inference_type = l.chomp.split("\t")
  tfs = cisbp_uniprot_acs.split(';').map{|uniprot_ac|
    request_uniprot_data_by_uniprot_ac(uniprot_ac)[:uniprot_id]
  }.uniq
  tfs.map{|tf|
    [tf, {cisbp_families: cisbp_families, cisbp_dbds: cisbp_dbds}]
  }
}.uniq.to_h

tf_to_cisbp.default_proc = ->(h,k){h[k] = {} }
# $stderr.puts(tf_to_cisbp)

jaspar_motifs = File.readlines('source_data/annotation/jaspar_prefinal_3.tsv').map{|l|
  motif, tf, gene, source_type, motif_length, family_1, family_2, family_3, family_4 = l.chomp.split("\t")
  ['jaspar', motif, tf, gene, source_type, motif_length, family_1, family_2, family_3, family_4, *tf_to_cisbp[tf].values_at(:cisbp_families, :cisbp_dbds)]
}

hocomoco_motifs = File.readlines('source_data/annotation/hocomoco_prefinal_2.tsv').map{|l|
  motif, tf, gene, source_type, motif_length, family_1, family_2, family_3, family_4 = l.chomp.split("\t")
  ['hocomoco', motif, tf, gene, source_type.sub('ChIP-Seq', 'ChIP-seq'), motif_length, family_1, family_2, family_3, family_4, *tf_to_cisbp[tf].values_at(:cisbp_families, :cisbp_dbds)]
}

cisbp_motifs = File.readlines('source_data/annotation/cisbp_refined.tsv').flat_map{|l|
  motif, motif_name, gene, cisbp_families, cisbp_dbds, cisbp_uniprot_acs, family_1, family_2, family_3, family_4, source_type, sr_model, inference_type = l.chomp.split("\t")
  motif_filename = (inference_type == 'direct') ? "ccg.epfl.ch/pwmtools/pwmlibs/cisbp/#{motif_name}.mat" : "ccg.epfl.ch/pwmtools/pwmlibs/cisbp_inf/#{motif_name}.mat"
  begin
    motif_length = File.readlines(motif_filename).map(&:strip).reject(&:empty?).reject{|l| l.start_with?('>') }.size
  rescue
    motif_length = nil
    $stderr.puts "Motif length for #{motif_name} not found"    
  end
  uniprot_acs = cisbp_uniprot_acs.split(';')
  if !uniprot_acs.empty?
    uniprot_acs.map{|uniprot_ac|
      uniprot_id = request_uniprot_data_by_uniprot_ac(uniprot_ac)[:uniprot_id]
      tf = uniprot_id && !uniprot_id.empty? ? uniprot_id : uniprot_ac
      ["cisbp_#{inference_type}", motif, tf, gene, source_type, motif_length, family_1, family_2, family_3, family_4, cisbp_families, cisbp_dbds]
    }.uniq
  else
    [["cisbp_#{inference_type}", motif, nil, gene, source_type, motif_length, family_1, family_2, family_3, family_4, cisbp_families, cisbp_dbds]]
  end
}

header = ['collection', 'motif', 'tf', 'gene', 'source_type', 'motif_length', *(1..4).map{|i| "TF_Class_level_#{i}"}, 'cisbp_families', 'cisbp_dbds']
puts header.join("\t")
[
  *jaspar_motifs,
  *hocomoco_motifs,
  *cisbp_motifs,
].each{|row| puts row.join("\t") }
