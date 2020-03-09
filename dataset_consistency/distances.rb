def norm(a)
  a.map{|ai| ai ** 2 }.sum(0.0) ** 0.5
end

def cosine(a, b)
  a.zip(b).map{|ai, bi| ai * bi }.sum(0.0) / (norm(a) * norm(b))
end

def normed_l1(a, b)
  a.zip(b).map{|ai, bi| (ai - bi).abs }.sum(0.0) / a.size
end

def normed_l2(a, b)
  (a.zip(b).map{|ai, bi| (ai - bi) ** 2 }.sum(0.0) / a.size) ** 0.5 
end


def l_inf(a, b)
  a.zip(b).map{|ai, bi| (ai - bi).abs }.max
end

def consistencies(matrix_1, matrix_2)
  datasets = (matrix_1.datasets & matrix_2.datasets)
  motifs = (matrix_1.motifs & matrix_2.motifs)

  matrix_2 = matrix_2.bound_at(datasets, motifs)
  matrix_1 = matrix_1.bound_at(datasets, motifs)

  datasets.map{|dataset|
    jy_vect = matrix_1.dataset_aucs(dataset)
    jolma_vect = matrix_2.dataset_aucs(dataset)
    info = {
      normed_l1: normed_l1(jy_vect, jolma_vect),
      normed_l2: normed_l2(jy_vect, jolma_vect),
      l_inf: l_inf(jy_vect, jolma_vect),
      cosine: cosine(jy_vect, jolma_vect)
    }.map{|k, v|
      [k, v.round(3)]
    }.to_h
    [dataset, info]
  }.to_h
end
