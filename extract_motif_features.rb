require 'bioinform'

class Array
  def mean; empty? ? nil : sum(0.0) / size; end
end

class Numeric
  def log_fact
    Math.lgamma(self + 1).first
  end
end

def dic(pos)
  n = pos.sum
  (pos.map{|el| el.log_fact }.sum - n.log_fact).to_f / n
end

def kdic(pos)
  n = pos.sum
  probs = pos.map{|el| el.to_f / n }
  background = [0.25] * 4
  dic(pos) - probs.zip(background).map{|prob, background_prob| prob * Math.log(background_prob) }.sum
end

def total_kdic(pcm)
  pcm.matrix.map{|pos| kdic(pos) }.sum
end

def information_content(pos)
  n = pos.sum
  probs = pos.map{|el| el.to_f / n }
  2 + probs.map{|prob| (prob == 0) ? 0 : prob * Math.log2(prob) }.sum
end

def mean_information_content(pcm)
  pcm.matrix.map{|pos| information_content(pos) }.mean
end

def gc_content(pos)
  (pos[1] + pos[2]).to_f / pos.sum
end

def total_gc_content(pcm)
  pcm.matrix.map{|pos| gc_content(pos) }.mean
end

filename = ARGV[0] #'source_data/annotation/motif_annotation_final.tsv'
header = File.readlines(filename).first.chomp.split("\t")
rows = File.readlines(filename).drop(1).map{|l|
  l.chomp.split("\t", header.size)
}

puts [*header, 'mean_IC', 'GC-content'].join("\t")

rows.drop(1).map{|row|
  motif = row[1]
  case row[0]
  when 'jaspar'
    pcm_fn = "source_data/motifs/jaspar_pcm/#{motif}.pcm"
  when 'hocomoco'
    pcm_fn = "source_data/motifs/hocomoco_pcm/#{motif}.pcm"
  when 'cisbp_direct'
    pcm_fn = Dir.glob("ccg.epfl.ch/pwmtools/pwmlibs/cisbp/#{motif}_*.mat").tap{|xs| raise if xs.size != 1 }.first
  when 'cisbp_inferred'
    pcm_fn = Dir.glob("ccg.epfl.ch/pwmtools/pwmlibs/cisbp_inf/#{motif}_*.mat").tap{|xs| raise if xs.size != 1 }.first
  end
  [row, pcm_fn]
}.each{|row, pcm_fn|
  # don't check counts discrepancy
  pcm = Bioinform::MotifModel::PCM.from_file(pcm_fn, validator: Bioinform::MotifModel::PM::VALIDATOR.make_strict)
  puts [*row, mean_information_content(pcm), total_gc_content(pcm)].join("\t")
}
