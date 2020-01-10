require 'json'
require 'httparty'
require 'WingenderTFClass'

require_relative 'gene2uniprot_curration'

HUMAN_RECOGNIZERS = WingenderTFClass::ProteinFamilyRecognizers::HumanAtLevel
MOUSE_RECOGNIZERS = WingenderTFClass::ProteinFamilyRecognizers::MouseAtLevel

def tf_classes_as_array(uniprot_id)
  (1..6).map{|level|
    [
      HUMAN_RECOGNIZERS[level].subfamilies_by_uniprot_id(uniprot_id),
      MOUSE_RECOGNIZERS[level].subfamilies_by_uniprot_id(uniprot_id),
    ].flatten
  }
end

TFCLASS_LEVELS = [:tf_superclass, :tf_class, :tf_family, :tf_subfamily, :tf_genus, :tf_molecular_species]

def tf_classes_as_hash(uniprot_id)
  TFCLASS_LEVELS.zip(tf_classes_as_array(uniprot_id)).to_h
end

def request_uniprot_data_by_gene_name(gene_name, organism:)
  if !$gene_cache
    if File.exist?('uniprot_cache_gene.json')
      $gene_cache = JSON.parse(File.read('uniprot_cache_gene.json')).transform_values{|v| v.map{|el| el.transform_keys(&:to_sym) } }
    else
      $gene_cache = {}
    end
  end
  return $gene_cache[gene_name]  if $gene_cache[gene_name]
  $gene_cache[gene_name] ||= begin
    base_url = "https://www.uniprot.org/uniprot/"
    params = {
      # organism:9606 - Homo sapiens (Human); Text query for organism obtains genes of other species 
      # like "Human adenovirus" organism on par with Human itself
      query: "gene_exact:#{gene_name}+organism:#{organism}+reviewed:yes",
      sort: 'score',
      columns: 'id,entry name,genes',
      # limit: 1,
      format: 'tab'
    }
    HTTParty.get(base_url, query: params).body.lines.drop(1).map{|l|
      uniprot_ac, uniprot_id, gene_names = l.chomp.split("\t", 3).map(&:strip)
      gene_names = gene_names.split(/\s+/).map(&:strip)
      {uniprot_ac: uniprot_ac, uniprot_id: uniprot_id, gene_names: gene_names}
    }
  end
  File.write('uniprot_cache_gene.json', $gene_cache.to_json)
  $gene_cache[gene_name]
end


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
    [tf, {cisbp_families: cisbp_families.split(';'), cisbp_dbds: cisbp_dbds.split(';')}]
  }
}.uniq.to_h
tf_to_cisbp.default_proc = ->(h,k){h[k] = {cisbp_families: [], cisbp_dbds: []} }


all_tf_names = [
  'source_data/motifs/hocomoco_genes2mat.txt',
  'source_data/motifs/jaspar_genes2mat.txt',
  'source_data/motifs/genes2cisbp.txt',
  'source_data/chipseq/remap_genes2exp.txt',
  'source_data/selex/jolma13_genes2exp.txt',
  'source_data/selex/genes2jolma_yang.txt',
  'source_data/uniprobe/genes2uniprobe.txt',
  'source_data/uniprobe/genes2uniprobe_manual.txt',
].flat_map{|fn|
  File.readlines(fn).map{|l|
    l.split("\t")[0].strip
  }
}.uniq.sort

###########################

$stderr.puts('Loading uniprot ids by TF names')

uniprot_tf_infos_unchecked = all_tf_names.reject{|tf_name|
  CURRATED_UNIPROT_AC.has_key?(tf_name) # we will add them later
}.map{|tf_name|
  $stderr.print('.')
  if tf_name.upcase == tf_name
    {tf_name: tf_name, uniprot_data: request_uniprot_data_by_gene_name(tf_name, organism: 9606)} # Homo sapiens
  else
    {tf_name: tf_name, uniprot_data: request_uniprot_data_by_gene_name(tf_name, organism: 10090)} # Mus musculus
  end
}

$stderr.puts("\nUniprots loaded")

###########################

# Detect motifs which can't be automatically assigned UniprotID to manually currate them (see gene2uniprot_curration.rb)
uniprot_tf_infos_unchecked.select{|info|
  info[:uniprot_data].size > 1
}.each{|info|
  info[:uniprot_data].each{|record|
    $stderr.puts [info[:tf_name], record].join("\t")
  }
}

uniprot_tf_infos_unchecked.select{|info|
  info[:uniprot_data].empty?
}.each{|info|
  $stderr.puts info[:tf_name]
}


######### Certain TF uniprots #########
uniprot_tf_infos_unchecked.select{|info|
  info[:uniprot_data].size != 1
}.each{|info|
  $stderr.puts "Probably incorrect mapping: #{info}"
}

uniprot_tf_infos = []
uniprot_tf_infos += uniprot_tf_infos_unchecked.select{|info|
  info[:uniprot_data].size == 1
}.map{|info|
  [info[:tf_name], info[:uniprot_data].first]
}
####### Uncertain TF uniprots (manual curation involved) #########
$stderr.puts('Loading uniprot ids by TF names from manually currated list')

uniprot_tf_infos += all_tf_names.select{|tf_name|
  CURRATED_UNIPROT_AC.has_key?(tf_name) # we will add them later
}.map{|tf_name|
  $stderr.print('.')
  uniprot_ac = CURRATED_UNIPROT_AC[tf_name]
  uniprot_record = request_uniprot_data_by_uniprot_ac(uniprot_ac)
  [tf_name, uniprot_record]
}

$stderr.puts("\nUniprots loaded")

#######################################
$stderr.puts("Load TFClass")
# uniprot + TFClass in a single hash
tf_infos = uniprot_tf_infos.map{|tf_name, uniprot_record|
  uniprot_id = uniprot_record[:uniprot_id]
  uniprot_ac = uniprot_record[:uniprot_ac]
  tfclass_hierarchy = tf_classes_as_hash(uniprot_id)
  cisbp_family_infos = tf_to_cisbp[uniprot_id]
  
  {
    tf_gene_name: tf_name,
    uniprot_id:uniprot_id, uniprot_ac: uniprot_ac,
    gene_names: uniprot_record[:gene_names],
    **tfclass_hierarchy,
    **cisbp_family_infos,
  }
}.sort_by{|info|
  info[:tf_name]
}

Homologene = Struct.new(:homologene_id, :organism_name, :taxon_id, :symbol, :entrezgene_id, :mgi_id, :hgnc_id, :omim_gene_id, :genetic_location, :genomic_coordinates, :nucl_refseq_ids, :protein_refseq_ids, :swissprot_ids) do
  def self.from_string(str)
    self.new(*str.chomp.split("\t"))
  end

  def self.each_in_file(filename, &block)
    return enum_for(:each_in_file, filename)  unless block_given?
    File.readlines(filename).map{|l|
      self.from_string(l)
    }.each(&block)
  end
end

human_symbols_by_mouse = Homologene.each_in_file('source_data/HOM_MouseHumanSequence.rpt').group_by(&:homologene_id).map{|homologene_id, records|
  mouse_records = records.select{|r| r.taxon_id == '10090' }
  human_records = records.select{|r| r.taxon_id == '9606' }
  human_symbols = human_records.map{|r| r.symbol }.uniq
  mouse_records.map{|r| [r.symbol, human_symbols] }
}.flatten(1).to_h

tf_infos_by_gene_name = tf_infos.group_by{|info| info[:tf_gene_name] }.to_h

tf_infos.select{|info|
  info[:uniprot_id].end_with?('MOUSE')
}.each{|info|
  mouse_gene_symbol = info[:tf_gene_name]
  human_gene_symbols = human_symbols_by_mouse[ mouse_gene_symbol ]
  human_infos = tf_infos_by_gene_name.values_at(*human_gene_symbols).flatten.compact
  human_cisbp_families = human_infos.map{|human_info| human_info[:cisbp_families]}.flatten
  human_cisbp_dbds = human_infos.map{|human_info| human_info[:cisbp_dbds]}.flatten
  info[:cisbp_families] = human_cisbp_families
  info[:cisbp_dbds] = human_cisbp_dbds
}; nil


$stderr.puts("TFClass loaded")

File.open('all_tf_infos.json', 'w') do |fw|
  tf_infos.each{|info|
    fw.puts(info.to_json)
  }
end


File.open('all_tf_infos.tsv', 'w') do |fw|
  header = ['tf_gene_name', 'uniprot_id', 'uniprot_ac', 'gene_names', TFCLASS_LEVELS, 'cisbp_families', 'cisbp_dbds']
  fw.puts(header.flatten.join("\t"))
  tf_infos.each{|info|
    row = info.values_at(:tf_gene_name, :uniprot_id, :uniprot_ac)
    row << info[:gene_names].join(':')
    row += info.values_at(*TFCLASS_LEVELS).map{|classes|
      classes.join(':')
    }
    row += info.values_at(:cisbp_families, :cisbp_dbds)
    fw.puts(row.join("\t"))
  }
end
