#!/usr/bin/env

import subprocess as sp
import gzip as gz
import optparse as opt
import cPickle
import vcf

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


def main():

    ## parse command line arguments
    parser = opt.OptionParser()
    parser.add_option('--vcf', dest = 'vcfFile', action = 'store', 
        help = 'File path for vcf for which to generate QC metrics.')
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
    parser.add_option('--rediscover-out', dest = 'kgOut', action = 'store',
        help = 'Output file location for rediscovery rate plots PDF.')
    parser.add_option('--hardy-out', dest = 'hardyOut', action = 'store',
        help = 'Output file location for Hardy Weinberg analysis plots PDF.')
    parser.add_option('--mendel-out', dest = 'mendelOut', action = 'store', 
        help = 'Output file location for Mendel inconsistency plots PDF.')
    parser.add_option('--temp-dir', dest = 'tempDir', action = 'store',
        help = 'Directory for writing intermediate analysis files.')
    (options, args) = parser.parse_args()

    vcfFile = options.vcfFile
    intermedBase = options.tempDir + 'qc'


    ## convert vcf to plink formatted files
    convertCall = ['vcftools', '--vcf', vcfFile, '--plink', '--recode', '--out',
        intermedBase]
    convertCall = smart_vcftools(convertCall)
    sp.call(convertCall)

    if options.pedFile:   
      ## pull pedigree information from previous run
      pedFile = options.vcfFile 
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
      plinkFile = intermedBase + '.pedigree.ped'
      plinkConn = open(plinkFile, 'w')
      with smart_open(intermedBase + '.ped') as plinkFile:
          for line in plinkFile:
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

    ## identify the total number of variants in the VCF file
    nsnp = len( [ rec for rec in vcf.Reader(open(vcfFile, 'r')) ] )

    ## get the ts/tv ratio
    tstvCall = ['vcftools', '--vcf', vcfFile, '--TsTv', nsnp, '--out',
        intermedBase + '.tstv']
    tstvCall = smart_vcftools(tstvCall)

    ## Calculate sample heterozygosity
    hetCall = ['vcftools', '--vcf', vcfFile, '--het', '--out', intermedBase ]
    hetCall = smart_vcftools(hetCall)
 
    ## calculate mendelian errors
    mendelCall = ['plink',
        '--ped', plinkFile,
        '--map', intermedBase + '.map',
        '--mendel',
        '--allow-no-sex',
        '--out', intermedBase]

    ## generate minor allele frequency distribution
    mafCall = ['plink',
    '--ped', intermedBase + '.pedigree.ped',
    '--map', intermedBase + '.map',
    '--freq',
    '--allow-no-sex',
    '--out', intermedBase]

    
    ## get the missingness rates
    missCall = ['plink',
    '--ped', intermedBase + '.pedigree.ped',
    '--map', intermedBase + '.map',
    '--missing',
    '--out', intermedBase]


    ## get hardy weinberg estimates
    hweCall = ['plink',
    '--ped', intermedBase + '.pedigree.ped',
    '--map', intermedBase + '.map',
    '--hardy',
    '--out', intermedBase]
    

    ## call external binaries and generate plots
    if options.rediscoverOut:
      ## calculate EVS rediscovery rate
      evsVar = cPickle.load(open('/nas40t0/vasya/exome_variant_server/ESP6500SI-V2-SSA137.snps.set.pickle', 'rb'))
      evsRes = calc_rediscovery(evsVar, vcfFile)
      del evsVar
      print >> open(intermedBase+'.evs', 'w'), '\t'.join(evsRes)

      ## calculate 1kG rediscovery rate
      kgVar = cPickle.load(open('/nas40t0/vasya/1kG/ALL_1000G_phase1integrated_v3_impute_var.pickle'))
      kgRes = calc_rediscovery(kgVar, vcfFile) 
      del kgVar
      print >> open(intermedBase+'.kg', 'w'), '\t'.join(kgRes)


    ## call external binaries and generate plots 
    if options.mendelOut:
      sp.call(mendelCall)
    if options.tstvOut:
      sp.call(tstvCall)
    if options.hetOut:
      sp.call(hetCall)
    if options.hardyOut:
      sp.call(hweCall)
    if options.missOut:
      sp.call(missCall)
    if options.mafOut:
      sp.call(mafCall)


if __name__ == '__main__':
    main()
