lines = ARGF.readlines.map{|l| l.chomp.split("\t") }
header = lines.first
row_header = lines.drop(1).map{|row| row.first }
matrix = lines.drop(1).map{|row| row.drop(1) }
columns = header.zip(matrix.transpose)
columns_sorted = columns.sort
header_sorted = columns_sorted.map(&:first)
matrix_sorted = columns_sorted.map(&:last).transpose
puts(header_sorted.join("\t"))
row_header.zip(matrix_sorted).each{|row_hdr, row|
  puts([row_hdr, *row].join("\t"))
}
