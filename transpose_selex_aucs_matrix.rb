filename = ARGV[0]
File.open(filename){|f|
  motifs = f.readline.chomp.split("\t")
  rows = f.readlines.map{|l| l.chomp.split("\t") }
  experiments = rows.map(&:first)
  matrix = rows.map{|row| row.drop(1) }
  puts experiments.join("\t")
  motifs.zip(matrix.transpose).each{|motif, row|
    puts [motif, *row].join("\t")
  }
}