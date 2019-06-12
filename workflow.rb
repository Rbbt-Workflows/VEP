#!/usr/bin/env ruby

require 'rbbt-util'
require 'rbbt/workflow'
require 'VEP'

Workflow.require_workflow "Sequence"

module VEP
  extend Workflow
  SOFTWARE_DIR=Rbbt.software.opt["ensembl-vep"].find
  CACHE_DIR=SOFTWARE_DIR.data

  #dep Sequence, :mutations_to_vcf
  #task :prepare => :array do |mutations|
  #  TSV.traverse step(:reference), :type => :array, :into => :stream do |line|
  #    next if line =~ /^#/

  #    mutation, ref = line.split "\t"
  #    next if ref.nil?
  #    chr, pos, mut = mutation.split(":")

  #    if %w(A C T G).include? mut
  #      start, eend = pos, pos
  #    elsif mut[0] == "-"

  #    end



  #    [chr, start, eend, ref+"/"+mut, '+']  * "\t"
  #  end
  #end

  extension :vcf
  input :args_VEP, :string, "Extra arguments for VEP"
  input :vcf_file, :text, "VCF File"
  input :organism, :string, "Organism code", nil
  input :reference, :string, "Reference code", nil
  dep Sequence, :mutations_to_vcf, :full_reference_sequence => true do |jobname,options|
    if options[:vcf_file]
      nil
    else
      {:input => options, :jobname => jobname}
    end
  end
  task :analysis => :text do |args_VEP,vcf_file,organism,reference|
    script = SOFTWARE_DIR["vep"].find

    reference ||=  Organism.hg_build organism
    grch = case reference
           when "hg19", 'b37'
             "GRCh37"
           when "hg38"
             "GRCh38"
           when "GRCh37", "GRCh38"
             reference
           else
             raise "Unknown build for #{organism} only GRCH37 GRCH38 allowed"
           end

    TmpFile.with_file do |tmpdir|
      Open.mkdir tmpdir
      vcf = TSV.get_stream(vcf_file || step(:mutations_to_vcf))
      CMD.cmd_log("perl #{script} --dir_plugins '/home/mvazque2/.vep/Plugins' --dir '#{CACHE_DIR}' --format vcf -o '#{tmpdir}/output' --quiet --assembly #{grch} --cache --offline --stats_text --force_overwrite --vcf --fork 20 #{args_VEP || ""}", :in => vcf)
      Open.read("#{tmpdir}/output")
    end
  end

  dep :analysis
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  task :mutated_isoforms => :tsv do |organism|
    dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Mutated Isoform"], :type => :flat, :namespace => organism
    dumper.init
    enst2enp = Organism.transcripts(organism).index :target => "Ensembl Protein ID", :fields => ["Ensembl Transcript ID"], :persist => true
    dumper = TSV.traverse step(:analysis), :type => :array, :into => dumper, :bar => "Processing VEP VCF" do |line|
      next if line =~ /^#/
      chr, pos, id, ref, alt, qual, filter, info = line.split("\t")
      next if info.nil?

      mutation = [chr, pos, alt] * ":"

      mis = info.split(',').collect do |str|
        next if str == '.'
        vals = str.split("|")
        enst = vals[6]
        next unless enst =~ /ENST/
        aa_pos = vals[14]
        next unless aa_pos =~ /\d/
        wt, alt = vals[15].split("/")
        alt = wt if alt.nil?
        
        ensp = enst2enp[enst]
        [ensp, [wt,aa_pos,alt]*""] * ":"
      end

      [mutation, mis]
    end
    io = Misc.collapse_stream(dumper.stream)
    CMD.cmd('tr "|" "\t"', :in => io, :pipe => true)
  end

  export_asynchronous :analysis, :mutated_isoforms
end
