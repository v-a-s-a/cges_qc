#!/usr/bin/python

from qc_pipeline import (make_plink_data,
                        make_mendel_data,
                        make_maf_data,
                        make_tstv_data,
                        make_miss_data,
                        compile_resource_descr,
                        get_base_dir)
import optparse as opt
import subprocess as sp
import os

def __main__():
  '''
  Generating comparisons across different levels of concordance.
  '''
  parser = opt.OptionParser()
  parser.add_option("--consensus-vcf", dest="consensus")
  parser.add_option("--3of4-vcf", dest="high")
  parser.add_option("--2of4-vcf", dest="mid")
  parser.add_option("--union-vcf", dest="union")
  parser.add_option("--ped-file", dest="pedFile")
  parser.add_option("--data-dir", dest="dataDir")
  parser.add_option("--mendel-out", dest="mendelOut")
  parser.add_option("--maf-out", dest="mafOut")
  parser.add_option("--miss-out", dest="missOut")
  parser.add_option("--tstv-out", dest="tstvOut")
  (options, args) = parser.parse_args()
  opt_dict = vars(options)
  fnull = open(os.devnull, 'w')

  ## store locations of resources
  branches = ['consensus', 'high', 'mid', 'union']
  resources = dict()
  for branch in branches:
    resources[branch] = compile_resource_descr(name = branch,
                                               vcf = opt_dict[branch],
                                               tmpdir = options.dataDir) 
  ## generate PLINK files if we are computing something that requires that information
  if options.mendelOut or options.mafOut:
    for branchData in resources.values():
      if not os.path.isfile(branchData['ped']):
        print "Generating PLINK files for: %s" % os.path.basename(branchData["vcf"])
        make_plink_data(vcfFile = branchData['vcf'],
                        pedFile = options.pedFile,
                        temp = branchData['temp'])
  
  if options.mendelOut:
    print "Generating Mendelian inconsistency plots."
    for branchData in resources.values():
      make_mendel_data(inputPed=branchData['ped'],
                       inputMap=branchData['map'],
                       temp=branchData['temp'])
    sp.call(['Rscript', get_base_dir() + '/R/concordance_mendel.R',
              resources['consensus']['temp'],
              resources['high']['temp'],
              resources['mid']['temp'],
              resources['union']['temp'],
              options.mendelOut])

  if options.mafOut:
    print "Generating minor allele frequency plots."
    for branchData in resources.values():
      make_maf_data(plinkMap=branchData['map'],
                    plinkPed=branchData['ped'],
                    temp=branchData['temp'])
    sp.call(['Rscript', get_base_dir() + '/R/concordance_maf.R',
            resources['consensus']['temp'],
            resources['high']['temp'],
            resources['mid']['temp'],
            resources['union']['temp'],
            options.mafOut])

  if options.tstvOut:
    print "Generating Ts/Tv plots."
    for branchData in resources.values():
      make_tstv_data(vcfFile=branchData['vcf'], temp=branchData['temp'])
    sp.call(['Rscript', get_base_dir() + '/R/concordance_tstv.R',
            resources['consensus']['temp'],
            resources['high']['temp'],
            resources['mid']['temp'],
            resources['union']['temp'],
            options.tstvOut], stdout = fnull)

  if options.missOut:
    print "Generating missingness plots."
    for branchData in resources.values():
      make_miss_data(plinkMap=branchData['map'],
                     plinkPed=branchData['ped'],
                     temp=branchData['temp'])
    sp.call(['Rscript', get_base_dir() + '/R/concordance_missing.R',
              resources['consensus']['temp'],
              resources['high']['temp'],
              resources['mid']['temp'],
              resources['union']['temp'],
              options.missOut])


if __name__ == "__main__":
  __main__()

