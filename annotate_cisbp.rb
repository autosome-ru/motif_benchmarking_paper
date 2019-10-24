require 'WingenderTFClass'


column_mapping = {
  motif_name: 'Motif_name', motif_id: 'Motif_ID', species: 'TF_Species', tf_id: 'TF_ID',
  family_id: 'Family_ID', tsource_id: 'TSource_ID', msource_id: 'MSource_ID', dbid: 'DBID',
  tf_name: 'TF_Name', tf_status: 'TF_Status', family_names: 'Family_Name', dbds: 'DBDs', dbd_count: 'DBD_Count', 
  cutoff: 'Cutoff', dbid_2: 'DBID.1',  motif_type: 'Motif_Type',
  msource_identifier: 'MSource_Identifier', msource_type: 'MSource_Type', pmids: 'PMID',
  msource_version: 'MSource_Version', sr_model: 'SR_Model', sr_no_threshold: 'SR_NoThreshold',
  tf_source_name: 'TfSource_Name', tfsource_url: 'TfSource_URL', 
  uniprot_swissprots: 'uniprot.Swiss-Prot', uniprot_tr_embls: 'uniprot.TrEMBL',
  inference_type: 'inference_type',
}

CisbpMotif = Struct.new('CisbpMotif', *column_mapping.keys, keyword_init: true) do
  COLUMN_MAPPING = column_mapping.freeze
  INVERT_COLUMN_MAPPING = column_mapping.invert.freeze
  def self.from_string(str, header, **kwargs)
    row = str.chomp.split("\t")
    info = header.map{|hdr| INVERT_COLUMN_MAPPING[hdr] }.zip(row).to_h
    info.merge!(kwargs)
    info[:sr_no_threshold] = (info[:sr_no_threshold] == 'None') ? nil : Integer(info[:sr_no_threshold])
    info[:motif_type] = (info[:motif_type] == 'None') ? nil : info[:motif_type]
    info[:dbd_count] = (info[:dbd_count] == 'None') ? nil : Integer(info[:dbd_count])
    info[:dbds] = info[:dbds].split(',')
    info[:family_names] = info[:family_names].split(',')
    info[:cutoff] = (info[:cutoff] == 'None') ? nil : Float(info[:cutoff])
    info[:pmids] = (info[:pmids] == 'None') ? [] : Array(info[:pmids])
    info[:msource_version] = (info[:msource_version] == 'None') ? nil : info[:msource_version]
    
    if info[:uniprot_swissprots] == 'None'
      info[:uniprot_swissprots] = []
    elsif info[:uniprot_swissprots].start_with?('[')
      info[:uniprot_swissprots] = info[:uniprot_swissprots][1...-1].split(',').map(&:strip).map{|x| x[1...-1] }
    else
      # do nothing
      info[:uniprot_swissprots] = Array(info[:uniprot_swissprots])
    end

    if info[:uniprot_tr_embls] == 'None'
      info[:uniprot_tr_embls] = []
    elsif info[:uniprot_tr_embls].start_with?('[')
      info[:uniprot_tr_embls] = info[:uniprot_tr_embls][1...-1].split(',').map(&:strip).map{|x| x[1...-1] }
    else
      info[:uniprot_tr_embls] = Array(info[:uniprot_tr_embls])
    end
    info[:inference_type] = info[:inference_type].to_sym
    self.new(**info)
  end

  def self.each_in_file(filename, **kwargs, &block)
    return enum_for(:each_in_file, filename, **kwargs)  unless block_given?
    File.open(filename){|f|
      header = f.readline.chomp.split("\t")
      f.each_line{|l|
        yield self.from_string(l, header, **kwargs)
      }
    }
  end
end

# We have only human-homologs
tf_classification_filename = WingenderTFClass::FilePaths::TFOntologyHuman
tf_classification = WingenderTFClass::OBO::TFClassification.from_file(tf_classification_filename)
RECOGNIZERS = Hash.new{|h, deepness|
  h[deepness] = WingenderTFClass::ProteinFamilyRecognizers::ByUniprotAC.new(tf_classification, deepness)
}

motifs = [
  *CisbpMotif.each_in_file('source_data/annotation/inferred_cisbp.tsv', inference_type: :inferred).to_a,
  *CisbpMotif.each_in_file('source_data/annotation/direct_cisbp.tsv', inference_type: :direct).to_a,
]
File.open('source_data/annotation/cisbp_refined.tsv', 'w'){|fw|
  motifs.each{|motif|
    uniprot_acs = motif.uniprot_swissprots
    families_levels = (1..4).map{|i| RECOGNIZERS[i].subfamilies_by_multiple_uniprot_acs(uniprot_acs) }
    
    # puts motif  if families.empty?
    info = [
      motif.motif_id, motif.motif_name, motif.tf_name, 
      motif.family_names.join(';'), motif.dbds.join(';'), motif.uniprot_swissprots.join(';'),
      *families_levels.map{|families| families.join(';') },
      motif.motif_type, motif.sr_model, motif.inference_type,
    ]
    fw.puts info.join("\t")
  }
}
