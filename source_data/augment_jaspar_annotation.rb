require 'json'
require 'httparty'
File.readlines('annotation/jasp_ann.tsv').each{|l|
  l.chomp!
  motif, gene, fam, subfam = l.split("\t", 4)
  matrix = motif.split('_').first
  begin
    opts = JSON.parse(HTTParty.get("http://jaspar.genereg.net/api/v1/matrix/#{matrix}/?format=json").body)
    species = opts['species'].map{|sp| sp['name'] }.join(';')
    type = opts['type']
    tf_class = opts['class'].compact.uniq.join(';')
    length = opts['pfm']['A'].size
    puts [motif, gene, fam, subfam, tf_class, species, type, length].join("\t")
  rescue
    $stderr.puts "Warning! Info for #{motif} can't be obtained"
    puts [motif, gene, fam, subfam, tf_class].join("\t")
  end
  # sleep(1)
}