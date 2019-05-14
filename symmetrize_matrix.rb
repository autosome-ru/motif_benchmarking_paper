header = readline
rows = readlines.map{|l|
  l.chomp.split("\t")
}

matrix_names = rows.map(&:first)
matrix = rows.map{|row| row.drop(1) }
resulting_matrix = Array.new(matrix.size){ Array.new(matrix.size) }
matrix.each_index{|i|
  matrix.each_index{|j|
    resulting_matrix[i][j] = (i >= j) ? matrix[i][j] : matrix[j][i]
  }
}

puts header
matrix_names.zip(resulting_matrix){|matrix_name, row|
  puts [matrix_name, *row].join("\t")
}
