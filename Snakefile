import itertools as it
import os

import numpy as np
import vcf
from pysradb.sraweb import SRAweb


SRRs = list(it.chain.from_iterable([
  [
    vcf_filename.split('.')[0]
    for vcf_filename
    in os.listdir(os.path.join('data', 'input', directory, 'Merged SnpEff vcf'))
  ] for directory in ['20200502', '20200509', '20200517']
]))


rule SRX:
  output:
    "data/srx.csv"
  run:
    db = SRAweb()
    db.srr_to_srx(SRRs).to_csv(output[0])

rule median_variant_distance:
  input:
    "data/vcf/{accession}.vcf"
  output:
    "data/vcf/{accession}-dist.txt"
  run:
    f = open(input[0], 'r')
    r = vcf.Reader(f)
    v = list(r)
    if len(v) > 0:
      p = [x.POS for x in v]
    else:
      p = [0]
    m = np.median(np.diff(p))
    with open(output[0], 'w') as o:
      o.write('%d' % m)
