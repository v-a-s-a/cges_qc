#!/bin/bash

## load virtual env
source ~/python2.7_env/bin/activate

python python/qc.pipeline.py \
  --atlas-vcf "test/atlas.test.vcf" \
  --gatk-vcf "test/gatk.test.vcf" \
  --freebayes-vcf "test/freebayes.test.vcf" \
  --cges-vcf "test/cges.test.vcf" \
  --ped-file "test/test.pedigree.txt" \
  --tstv-out "test/tstv.dat" \
  --het-out "test/het.dat" \
  --maf-out "test/maf.dat" \
  --miss-out "test/miss.dat" \
  --mendel-out "test/mendel.dat" \
  --temp-dir "test/" \
  --rediscover-out tests/rediscover.dat 
