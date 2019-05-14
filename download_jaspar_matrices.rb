require 'fileutils'
require 'httparty'
FileUtils.mkdir_p 'source_data/motifs/jaspar_pcm'
matrix_names = File.readlines('source_data/motifs/jaspar_genes2mat.txt').map{|l| l.chomp.split("\t").last }
matrix_names.shuffle.each{|matrix|
  begin
    output_fn = "source_data/motifs/jaspar_pcm/#{matrix}.pcm"
    next  if File.exist?(output_fn)

    matrix_id = matrix.split('_').first
    pcm = HTTParty.get("http://jaspar.genereg.net/api/v1/matrix/#{matrix_id}.pfm").body.lines
    mat = pcm.drop(1).map{|l|
      l.strip.split(/\s+/)
    }.transpose
    File.open(output_fn, 'w'){|fw|
      fw.puts ">#{matrix}"
      mat.each{|row|
        fw.puts row.join("\t")
      }
    }
    $stderr.print('.')
    sleep(0.2)
  rescue
    $stderr.puts "\n#{matrix} failed to be downloaded"
  end
}
