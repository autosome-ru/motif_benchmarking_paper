filename = ARGV[0]
normalize_cisbp = ARGV.delete('--normalize-cisbp')

File.open(filename){|f|
  motifs = f.readline.chomp.split("\t")
  motifs = motifs.map{|nm| nm.split('_',3)[0,2].join('_') } if normalize_cisbp
  rows = f.readlines.map{|l| l.chomp.split("\t") }
  experiments = rows.map(&:first)
  matrix = rows.map{|row| row.drop(1) }
  puts experiments.join("\t")
  motifs.zip(matrix.transpose).each{|motif, row|
    puts [motif, *row].join("\t")
  }
}
