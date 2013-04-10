#!/usr/bin/python
# Parses InParanoid results. This doesn't include converting
# gene names.
# This would be quite doable in Perl, or probably even R,
# but I wanted to try Python.


in_paranoid_dir = '/home/jburdick/gcb/data/homology/InParanoid/'

import gzip
import xml.etree.ElementTree as ET

outfile = gzip.open('/home/jburdick/gcb/git/data/homology/InParanoid.tsv.gz', 'w')
outfile.write('file\tcluster\tspecies\tgene\tscore\n')

# Parses one InParanoid XML file, and writes out
# some of the information in tab-separated format.
def parse_file(f) :

  tree = ET.parse(in_paranoid_dir + f)
  cl = tree.getroot().find('CLUSTERS')

  for cluster in cl:
    clusterno = cluster.attrib.get('CLUSTERNO')
    bitscore = cluster.attrib.get('BITSCORE')
#  print(clusterno + "  " + bitscore)
    for gene in cluster:
      geneid = gene.attrib.get('GENEID')
      score = gene.attrib.get('SCORE')
      species = gene.attrib.get('SPECIES')
      outfile.write(f + '\t' + clusterno + '\t' +
        species + '\t' + geneid + '\t' + score + '\n')

parse_file('InParanoid.C.elegans-H.sapiens.xml')
parse_file('InParanoid.C.elegans-M.musculus.xml')
parse_file('InParanoid.C.elegans-D.melanogaster.xml')

outfile.close()


