#!/usr/bin/env

import subprocess as sp
import gzip as gz
import optparse as opt
import cPickle
import vcf
import os

def smart_open(f):
    '''
    Open compressed files if compressed
    '''
    if f.endswith('.gz'):
        return gz.open(f)
    else:
        return open(f)

def smart_vcftools(call):
    '''
    Analog to smart open for calling vcftools with a gzipped vcf
    '''
    for idx, element in enumerate(call):
        if element.endswith('.vcf.gz'):
            call[idx-1] = '--gzvcf'
    return call

def makeID(rec, altall):
  var = rec.CHROM.lstrip('chr') + ':' + str(rec.POS) + '.' + rec.REF + '.' + str(altall)
  altvar = rec.CHROM.lstrip('chr') + ':' + str(rec.POS) + '.' + str(altall) + '.' + rec.REF
  return var, altvar

def calc_rediscovery(evs, vcfFile):
  '''
  Calculate the variant rediscovery rate in a VCF file.
  '''
  count = 0
  total = 0
  vcfReader = vcf.Reader(open(vcfFile))
  for rec in vcfReader:
    for allele in rec.ALT:
      var, altvar = makeID(rec, allele)
      if (var in evs) or (altvar in evs):
        count += 1
      total += 1
  return (vcfFile.split('/')[-1], str(float(count)/float(total)))

def make_plink_data(vcfFile, pedFile, temp, pedOut=None, mapOut=None):
  '''
  Given a VCF file, recode into PLINK files and cleanup auxilliary files.

  If a pedigree file is supplied, the recode the ped file with that information.

  '''
  fnull = open(os.devnull, 'w')

  ## convert vcf to plink formatted files
  convertCall = ['vcftools', '--vcf', vcfFile, '--plink', '--recode', '--out',
      temp ]
  convertCall = smart_vcftools(convertCall)
  sp.call(convertCall, stdout=fnull)

  if pedFile:   
    ## pull pedigree information from previous run
    mothers = {}
    fathers = {}
    families = {}
    sexes = {}
    with open(pedFile) as fconn:
      for line in fconn:
        line = line.split(',')
        mother = line[3]
        father = line[2]
        if mother=='Not Sequenced':
            mother = '0'
        if father=='Not Sequenced':
            father = '0'
        fathers[line[1]] = father
        mothers[line[1]] = mother
        family = line[0] 
        families[line[1]] = family
        sexes[line[1]] = line[4]

    ## add pedigree information to plink ped file
    plinkFile = temp + '.pedigree.ped'
    plinkConn = open(plinkFile, 'w')
    with smart_open(temp + '.ped') as pedPlinkFile:
        for line in pedPlinkFile:
            line = line.strip().split('\t')
            iid = line[1]
            if fathers.get(iid):
                line[2] = fathers[iid]
            if mothers.get(iid):
                line[3] = mothers[iid]
            if families.get(iid):
                line[0] = families[iid]
            if sexes.get(iid):
                line[4] = sexes[iid]
            ## add a fake phenotype
            line[5] = '2'
            print >> plinkConn, '\t'.join(line)
    plinkConn.close()

    ## move outputs and clean up vestigial files
    sp.call(['rm', temp+'.recode.vcf'])
    sp.call(['rm', temp+'.ped'])
    sp.call(['mv', temp+'.pedigree.ped', temp+'.ped'])
    if pedOut and mapOut:
      sp.call(['mv', temp+'.pedigree.ped', pedout])
      sp.call(['mv',  temp+'.map', mapOut])
    


def make_mendel_data(inputPed, inputMap, temp, trioOut=None, locusOut=None):
  '''
  Given a set of PLINK files, generate mendel stats, and move files to specified location.
  '''
  fnull = open(os.devnull, 'w')
  ## calculate mendelian errors
  mendelCall = ['plink',
    '--ped', inputPed,
    '--map', inputMap,
    '--mendel',
    '--allow-no-sex',
    '--out', temp]
  sp.call(mendelCall, stdout = fnull)

  ## clean up and move files
  sp.call(['rm', temp+'.log'])
  sp.call(['rm', temp+'.imendel'])
  sp.call(['rm', temp+'.nosex'])
  if trioOut and locusOut:
    sp.call(['mv', temp+'.fmendel', trioOut])
    sp.call(['mv', temp+'.lmendel', locusOut])


def make_tstv_data(vcfFile, temp, tstvOut=None):
  fnull = open(os.devnull, 'w')
  ## identify the total number of variants in the VCF file
  nsnp = len( [ rec for rec in vcf.Reader(open(vcfFile, 'r')) ] )

  ## get the ts/tv ratio
  tstvCall = ['vcftools', '--vcf', vcfFile, '--TsTv', str(nsnp), '--out',
      temp]
  tstvCall = smart_vcftools(tstvCall)
  sp.call(tstvCall, stdout = fnull) 
 
  ## clean up files
  sp.call(['rm', temp+'.log'])
  sp.call(['rm', temp+'.TsTv'])
  if tstvOut:
    sp.call(['mv', temp+'.TsTv.summary', tstvOut])

def make_het_data(vcfFile, temp, hetOut=None):
  fnull = open(os.devnull, 'w')
  ## Calculate sample heterozygosity
  hetCall = ['vcftools', '--vcf', vcfFile, '--het', '--out', temp ]
  hetCall = smart_vcftools(hetCall)
  sp.call(hetCall, stdout = fnull)

  ## clean up files
  sp.call(['rm', temp+'.log'])
  if hetOut:
    sp.call(['mv', temp+'.het', hetOut])

def make_maf_data(plinkMap, plinkPed, temp, mafOut=None):
  fnull = open(os.devnull, 'w')
  ## generate minor allele frequency distribution
  mafCall = ['plink',
    '--ped', plinkPed,
    '--map', plinkMap,
    '--freq',
    '--allow-no-sex',
    '--out', temp]
  sp.call(mafCall, stdout=fnull)

  ## clean up files
  sp.call(['rm', temp+'.log'])
  sp.call(['rm', temp+'.nosex'])
  if mafOut:
    sp.call(['mv', temp+'.frq', mafOut])


def make_miss_data(plinkMap, plinkPed, temp, missOut=None): 
  fnull = open(os.devnull, 'w')
  ## get the missingness rates
  missCall = ['plink',
    '--ped', plinkPed,
    '--map', plinkMap,
    '--missing',
    '--out', temp]
  sp.call(missCall, stdout = fnull) 
 
  ## clean up files
  sp.call(['rm', temp+'.log'])
  sp.call(['rm', temp+'.nosex'])
  if missOut:
    sp.call(['mv', temp+'.lmiss', missOut])


def make_hardy_data(plinkMap, plinkPed, temp, hardyOut=None):
  fnull = open(os.devnull, 'w')
  ## get hardy weinberg estimates
  hweCall = ['plink',
    '--ped', plinkPed,
    '--map', plinkMap,
    '--hardy',
    '--out', temp]
  sp.call(hweCall, stdout = fnull)

  ## clean up files
  sp.call(['rm', temp+'.log'])
  sp.call(['rm', temp+'.nosex'])
  if hardyOut:
    sp.call(['mv', temp+'.hwe', hardyOut])
    

def make_rediscovery_data(vcfFile, evsOut, kgOut):
  fnull = open(os.devnull, 'w')
  ## calculate EVS rediscovery rate
  evsVar = cPickle.load(open('/nas40t0/vasya/exome_variant_server/ESP6500SI-V2-SSA137.snps.set.pickle', 'rb'))
  evsRes = calc_rediscovery(evsVar, vcfFile)
  del evsVar
  print >> open(evsOut, 'w'), '\t'.join(evsRes)

  ## calculate 1kG rediscovery rate
  kgVar = cPickle.load(open('/nas40t0/vasya/1kG/ALL_1000G_phase1integrated_v3_impute_var.pickle'))
  kgRes = calc_rediscovery(kgVar, vcfFile) 
  del kgVar
  print >> open(kgOut, 'w'), '\t'.join(kgRes)


def compile_resource_descr(name, vcf, tmpdir):
  '''
  Make a dictionary of the file descriptions associated with a VCF file.
  '''
  resources = { "vcf": vcf,
                "ped": tmpdir + name + '.ped' ,
                "map": tmpdir + name + '.map',
                "temp": tmpdir + name }
  return resources

def main():
  ## parse command line arguments
  parser = opt.OptionParser()
  parser.add_option('--cges-vcf', dest = 'cges', action = 'store', 
      help = 'File path for CGES VCF for which to generate QC metrics.')
  parser.add_option('--atlas-vcf', dest = 'atlas', action = 'store', 
      help = 'File path for ATLAS VCF for which to generate QC metrics.')
  parser.add_option('--gatk-vcf', dest = 'gatk', action = 'store', 
      help = 'File path for GATK VCF for which to generate QC metrics.')
  parser.add_option('--freebayes-vcf', dest = 'freebayes', action = 'store', 
      help = 'File path for Freebayes VCF for which to generate QC metrics.')
  parser.add_option('--ped-file', dest = 'pedFile', action = 'store', 
      help = 'Pedigree file for samples (Optional).')
  parser.add_option('--tstv-out', dest = 'tstvOut', action = 'store', 
      help = 'Output file location for TsTv plots PDF.')
  parser.add_option('--het-out', dest = 'hetOut', action = 'store', 
      help = 'Output file location for heterozygosity plots PDF.')
  parser.add_option('--maf-out', dest = 'mafOut', action = 'store', 
      help = 'Output file location for minor allele frequency plots PDF.')
  parser.add_option('--miss-out', dest = 'missOut', action = 'store', 
      help = 'Output file location for missingess plots PDF.')
  parser.add_option('--rediscover-out', dest = 'rediscoverOut', action = 'store',
      help = 'Output file location for rediscovery rate plots PDF.')
  parser.add_option('--mendel-out', dest = 'mendelOut', action = 'store', 
      help = 'Output file location for Mendel inconsistency plots PDF.')
  parser.add_option('--temp-dir', dest = 'tempDir', action = 'store',
      help = 'Directory for writing intermediate analysis files.')
  (options, args) = parser.parse_args()
  opt_dict = vars(options)
  
  fnull = open(os.devnull, 'w')

  ## give each branch a name
  branches = ['atlas', 'gatk', 'freebayes', 'cges']

  ## store resource locations for each branch
  resources = dict()
  for branch in branches:
    resources[branch] = compile_resource_descr(name = branch, vcf = opt_dict[branch], tmpdir = options.tempDir)

  ## check if we need to recode into PLINK files
  if options.mendelOut or options.mafOut or options.missOut or options.hardyOut:
    for branchData in resources.values():
      make_plink_data(vcfFile=branchData['vcf'], pedFile=options.pedFile, temp=branchData['temp'])

  if options.mendelOut:
    ## generate data
    for branchData in resources.values():
      make_mendel_data(inputPed=branchData['ped'], inputMap=branchData['map'], temp=branchData['temp'])
    ## make plots PDF by executing an R script with passed in arguments
    sp.call(['Rscript', 'R/mendel.R', resources['atlas']['temp'], resources['gatk']['temp'], resources['freebayes']['temp'], resources['cges']['temp'], options.mendelOut], stdout = fnull)

  if options.tstvOut:
    ## generate data
    for branchData in resources.values():
      make_tstv_data(vcfFile=branchData['vcf'], temp=branchData['temp'])
    ## generate plot PDF
    sp.call(['Rscript', 'R/tstv.R', resources['atlas']['temp'], resources['gatk']['temp'], resources['freebayes']['temp'], resources['cges']['temp'], options.tstvOut], stdout = fnull)

  if options.hetOut:
    for branchData in resources.values():
      make_het_data(vcfFile=branchData['vcf'], temp=branchData['temp'])
    sp.call(['Rscript', 'R/het.R', resources['atlas']['temp'], resources['gatk']['temp'], resources['freebayes']['temp'], resources['cges']['temp'], options.hetOut], stdout = fnull)

  if options.mafOut:
    for branchData in resources.values():
      make_maf_data(plinkMap=branchData['map'], plinkPed=branchData['ped'], temp=branchData['temp'])
    sp.call(['Rscript', 'R/maf.R', resources['atlas']['temp'], resources['gatk']['temp'], resources['freebayes']['temp'], resources['cges']['temp'], options.mafOut])
     
  if options.missOut:
    for branchData in resources.values(): 
      make_miss_data(plinkMap=branchData['map'], plinkPed=branchData['ped'], temp=branchData['temp'])
    sp.call(['Rscript', 'R/missing.R', resources['atlas']['temp'], resources['gatk']['temp'], resources['freebayes']['temp'], resources['cges']['temp'], options.missOut])

#  if options.rediscoverOut:
#    make_rediscovery_data()
    



if __name__ == '__main__':
    main()
