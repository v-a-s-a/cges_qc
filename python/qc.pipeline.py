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
            print element
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
    parser.add_option('--intermed-base', dest = 'intermedBase', action = 'store', 
        help = 'File path base for storage of intermediate files.')
    (options, args) = parser.parse_args()

    vcfFile = options.vcfFile
    intermedBase = options.intermedBase

    ## identify the total number of variants in the VCF file
    zgrepCall = ['zgrep', '-v', '\'#\'', vcfFile]
    wcCall = ['wc', '-l']
    zgrep = sp.Popen( zgrepCall, stdout = sp.PIPE )
    wc = sp.Popen( wcCall, stdin = zgrep.stdout, stdout = sp.PIPE )
    zgrep.stdout.close()
    nsnp = wc.communicate()[0]

    ## get the ts/tv ratio
    tstvCall = ['vcftools', '--vcf', vcfFile, '--TsTv', nsnp, '--out',
        intermedBase + '.tstv']
    tstvCall = smart_vcftools(tstvCall)

    ## Calculate sample heterozygosity
    hetCall = ['vcftools', '--vcf', vcfFile, '--het', '--out', intermedBase ]
    hetCall = smart_vcftools(hetCall)

    ## convert vcf to plink formatted files
    convertCall = ['vcftools', '--vcf', vcfFile, '--plink', '--recode', '--out',
        intermedBase]
    convertCall = smart_vcftools(convertCall)
    sp.call(convertCall)
   
    ## pull pedigree information from previous run
    pedFile = '/nas40t0/vasya/autism/data/ha_vcf_trios_pedfile_012813_nick_vcf_trios.csv'
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
    newPlinkFile = intermedBase + '.pedigree.ped'
    newPlinkConn = open(newPlinkFile, 'w')
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
            print >> newPlinkConn, '\t'.join(line)
    newPlinkConn.close()
 
    ## calculate mendelian errors
    mendelCall = ['plink',
        '--ped', newPlinkFile,
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

    ## call plink and vcftools based statistics
    sp.call(mendelCall)
    sp.call(tstvCall)
    sp.call(hetCall)
    sp.call(hweCall)
    sp.call(missCall)
    sp.call(mafCall)


if __name__ == '__main__':
    main()
