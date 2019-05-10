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

def request_uniprot_data_by_gene_name(gene_name)
  base_url = "https://www.uniprot.org/uniprot/"
  params = {
    # organism:9606 - Homo sapiens (Human); Text query for organism obtains genes of other species 
    # like "Human adenovirus" organism on par with Human itself
    query: "gene_exact:#{gene_name}+organism:9606+reviewed:yes",
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

def request_uniprot_data_by_uniprot_ac(uniprot_ac)
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

all_tf_names = [
  'source_data/motifs/hocomoco_genes2mat.txt',
  'source_data/motifs/jaspar_genes2mat.txt',
  'source_data/chipseq/remap_genes2exp.txt',
  'source_data/selex/jolma13_genes2exp.txt',
].flat_map{|fn|
  File.readlines(fn).map{|l|
    l.split("\t")[0].strip
  }
}.sort.uniq

###########################

$stderr.puts('Loading uniprot ids by TF names')

uniprot_tf_infos_unchecked = all_tf_names.reject{|tf_name|
  CURRATED_UNIPROT_AC.has_key?(tf_name) # we will add them later
}.map{|tf_name|
  $stderr.print('.')
  {tf_name: tf_name, uniprot_data: request_uniprot_data_by_gene_name(tf_name)}
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
  
  {
    tf_gene_name: tf_name,
    uniprot_id:uniprot_id, uniprot_ac: uniprot_ac,
    gene_names: uniprot_record[:gene_names],
    **tfclass_hierarchy
  }
}.sort_by{|info|
  info[:tf_name]
}
$stderr.puts("TFClass loaded")

File.open('all_tf_infos.json', 'w') do |fw|
  tf_infos.each{|info|
    fw.puts(info.to_json)
  }
end


File.open('all_tf_infos.tsv', 'w') do |fw|
  header = ['tf_gene_name', 'uniprot_id', 'uniprot_ac', 'gene_names', TFCLASS_LEVELS]
  fw.puts(header.flatten.join("\t"))
  tf_infos.each{|info|
    row = info.values_at(:tf_gene_name, :uniprot_id, :uniprot_ac)
    row << info[:gene_names].join(':')
    row += info.values_at(*TFCLASS_LEVELS).map{|classes|
      classes.join(':')
    }
    fw.puts(row.join("\t"))
  }
end
