require 'WingenderTFClass'
recs = WingenderTFClass::ProteinFamilyRecognizers::HumanAtLevel
field = Integer(ARGV[0]) - 1
$stdin.each_line{|l|
  row = l.chomp.split("\t")
  uniprot_id = row[field]
  tf_superclass = recs[1].subfamilies_by_uniprot_id(uniprot_id)
  tf_class = recs[2].subfamilies_by_uniprot_id(uniprot_id)
  fams = recs[3].subfamilies_by_uniprot_id(uniprot_id)
  subfams = recs[4].subfamilies_by_uniprot_id(uniprot_id)

  classification = [tf_superclass, tf_class, fams, subfams].map{|cl_values| cl_values.compact.uniq.join(';') }
  infos = [*row, *classification]
  puts infos.join("\t")
}
