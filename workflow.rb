#!/usr/bin/env ruby

require 'rbbt-util'
require 'rbbt/workflow'
require 'VEP'

Workflow.require_workflow "Sequence"

module VEP
  extend Workflow
  SOFTWARE_DIR=Rbbt.software.opt["ensembl-tools"].produce["scripts/variant_effect_predictor"].find

  dep Sequence, :reference
  task :prepare => :array do |mutations|
    TSV.traverse step(:reference).grace, :type => :array, :into => :stream do |line|
      next if line =~ /^#/
        mutation, ref = line.split "\t"
      chr, pos, mut = mutation.split(":")
      [chr, pos, pos, ref+"/"+mut, '+']  * "\t"
    end
  end

  dep :prepare
  task :analysis => :tsv do 
    script = SOFTWARE_DIR["variant_effect_predictor.pl"].find
    data_dir = Rbbt.software.opt["ensembl-tools"]["Data.GRCh37"].find
    log :VEP, "Running VEP script"
    CMD.cmd("perl #{script} -i '#{step(:prepare).join.path}' -o '#{path}' --dir '#{data_dir}' --cache --offline --stats_text")
    nil
  end

  dep :analysis
  task :mutated_isoforms => :tsv do 
    organism = Organism.default_code("Hsa")
    dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Mutated Isoform"], :type => :flat, :organism => organism
    dumper.init
    enst2enp = Organism.transcripts(organism).index :target => "Ensembl Protein ID", :fields => ["Ensembl Transcript ID"], :persist => true
    dumper = TSV.traverse step(:analysis), :type => :array, :into => dumper do |line|
      next if line =~ /^#/
      mut_str, ehr_pos, ref, ensg, enst, type, conseq, cdna_pos, cds_pos, pos, change_str, *rest = line.split "\t", -1
      next if change_str.nil? or change_str.empty? or change_str == '-'
      mutation = mut_str.gsub("_",':').gsub(/:[^:]*?\//, ':')
      protein = enst2enp[enst]

      if change_str.length == 1
        change = [change_str, change_str] * pos.to_s
      elsif change_str.length > 3
        if conseq =~ /frameshift/
          change = [change_str.split("/").first, "FrameShift"] * pos.to_s
        else
          change = [change_str.split("/").first, "Indel"] * pos.to_s
        end
      else
        change = change_str.split("/") * pos.to_s
      end

      mi = protein + ':' + change
      [mutation, mi]
    end
    io = Misc.collapse_stream(dumper.stream)
    CMD.cmd('tr "|" "\t"', :in => io, :pipe => true)
  end
end
