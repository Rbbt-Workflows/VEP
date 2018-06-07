require 'rbbt-util'
require 'rbbt/resource'

module VEP
  def self.organism(org="Hsa")
    Organism.default_code(org)
  end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  Rbbt.claim Rbbt.software.opt["ensembl-tools"], :install, Rbbt.share.install.software.ENSEMBL.find(:lib)
  Rbbt.claim Rbbt.software.opt["ensembl-vep"], :install, Rbbt.share.install.software.VEP.find(:lib)
end
