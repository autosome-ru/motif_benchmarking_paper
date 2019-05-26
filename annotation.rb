require 'json'

TFCLASS_LEVELS = [:tf_superclass, :tf_class, :tf_family, :tf_subfamily, :tf_genus, :tf_molecular_species]

# "Some TF family{1.2.3}" --> [1,2,3]
def tfclass_id(tfclass_name)
  m = tfclass_name.match(/\{([\d.]+)\}/)
  m && m[1].split('.').map(&:to_i)
end

class Annotation
  attr_reader :experiments, :motifs
  # Specify subsets of experiments and motifs
  def initialize(experiments, motifs)
    @experiments = experiments
    @motifs = motifs
  end

  ############################

  def motif_tfs
    @motif_tfs ||= motifs.flat_map{|motif|
      tfs_by_motif(motif)
    }.uniq
  end

  def experiment_tfs
    @experiment_tfs ||= experiments.map{|experiment|
      tf_by_experiment(experiment)
    }.compact.uniq
  end

  ############################

  def motif_source_type(motif)
    @motif_type_map ||= begin
      motif_type_pairs = []
      motif_type_pairs += File.readlines('source_data/motifs/jaspar_infos.json').map{|l|
        JSON.parse(l.chomp)
      }.map{|info|
        ["#{info['matrix_id']}_#{info['name']}", info['type']]
      }
      motif_type_pairs += File.readlines('source_data/annotation/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv').drop(1).map{|l|
        row = l.chomp.split("\t")
        [row[0], row[7]]
      }

      motif_type_pairs  = motif_type_pairs.map{|k,v|
        [k, v&.sub('ChIP-Seq', 'ChIP-seq')]
      }
      motif_type_pairs.to_h
    end
    @motif_type_map[motif]
  end

  ############################

  def tf_motif_pairs
    @tf_motif_pairs ||= [
      'source_data/motifs/hocomoco_genes2mat.txt',
      'source_data/motifs/jaspar_genes2mat.txt',
    ].flat_map{|fn|
      File.readlines(fn).map{|l|
        l.split("\t").map(&:strip)
      }
    }.select{|tf, motif|
      motifs.include?(motif)
    }
  end

  def tfs_by_motif(motif)
    @tfs_by_motif ||= tf_motif_pairs.group_by{|tf, motif|
      motif
    }.map{|motif, pairs|
      [motif, pairs.map{|tf, motif| tf }]
    }.to_h
    @tfs_by_motif[motif] || []
  end

  def motifs_by_tf(tf)
    @motifs_by_tf ||= tf_motif_pairs.group_by{|tf, motif|
      tf
    }.map{|tf, pairs|
      [tf, pairs.map{|tf, motif| motif }]
    }.to_h
    @motifs_by_tf[tf] || []
  end

  ############################

  def experiment_tf_pairs
    @experiment_tf_pairs ||= [
      'source_data/chipseq/remap_genes2exp.txt',
      'source_data/selex/jolma13_genes2exp.txt',
    ].flat_map{|experiment_mapping_fn|
      File.readlines(experiment_mapping_fn).map{|l|
        tf, experiment = l.chomp.split("\t").first(2)
        [experiment, tf]
      }
    }
  end

  def tf_by_experiment(experiment)
    @tf_by_experiment ||= experiment_tf_pairs.to_h
    @tf_by_experiment[experiment]
  end

  def experiments_by_tf(tf)
    @experiments_by_tf ||= experiment_tf_pairs.select{|experiment, tf|
      # include only experiments in specified subset
      experiments.include?(experiment)
    }.group_by{|experiment, tf|
      tf
    }.map{|tf, exp_tf_pairs|
      [tf, exp_tf_pairs.map{|exp,tf| exp }]
    }.to_h
    @experiments_by_tf[tf] || []
  end

  #############################
  
  def tfclass_infos
    @tfclass_infos ||= File.readlines('all_tf_infos.json').map{|l|
      JSON.parse(l, symbolize_names: true)
    }
  end

  def tfclass_names(tfclass_level)
    @tfclass_names ||= {}
    @tfclass_names[tfclass_level] ||= tfclass_infos.flat_map{|info|
      info[tfclass_level]
    }.uniq.sort_by{|fam|
      tfclass_id(fam)
    }
  end

  def tfclass_parent(tfclass_name)
    child_id = tfclass_id(tfclass_name)
    parent_id = child_id[0...-1]
    parent_level = TFCLASS_LEVELS[ parent_id.length - 1 ]
    tfclass_names(parent_level).detect{|candidate_tfclass_name|
      tfclass_id(candidate_tfclass_name) == parent_id
    }
  end

  # def tfclass_families; tfclass_names(:tf_family); end
  # def tfclass_classes; tfclass_names(:tf_class); end

  def tf_infos_by_tfclass_name(tfclass_level, tfclass_name)
    @tf_infos_by_tfclass_name ||= {}
    @tf_infos_by_tfclass_name[tfclass_level] ||= begin
      hsh = Hash.new{|h,k| h[k] = [] }
      tfclass_infos.each{|info|
        info[tfclass_level].each{|tfclass_name|
          hsh[tfclass_name] << info
        }
      }
      hsh
    end
    @tf_infos_by_tfclass_name[tfclass_level][tfclass_name] || []
  end

  def tfs_by_tfclass_name(tfclass_level, tfclass_name)
    tf_infos_by_tfclass_name(tfclass_level, tfclass_name).map{|info|
      info[:tf_gene_name]
    }.compact.uniq
  end

  def motifs_by_tfclass_name(tfclass_level, tfclass_name)
    tfs_by_tfclass_name(tfclass_level, tfclass_name).flat_map{|tf|
      motifs_by_tf(tf)
    }.compact.uniq
  end

  #############################

  def tf_info_by_gene_name(name)
    @tf_info_by_gene_name ||= tfclass_infos.group_by{|info|
      info[:tf_gene_name]
    }.map{|k,vs|
      [k, vs.first]
    }.to_h
    @tf_info_by_gene_name[name]
  end

  #############################
  
  def experiments_by_tfclass_name(tfclass_level, tfclass_name)
    @experiments_by_tfclass_name ||= {}
    @experiments_by_tfclass_name[tfclass_level] ||= {}
    @experiments_by_tfclass_name[tfclass_level][tfclass_name] ||= begin
      tfs_by_tfclass_name(tfclass_level, tfclass_name).flat_map{|tf|
        experiments_by_tf(tf)
      }.uniq
    end
    @experiments_by_tfclass_name[tfclass_level][tfclass_name]
  end
end
