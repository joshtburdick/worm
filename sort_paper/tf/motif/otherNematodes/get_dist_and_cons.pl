#!/usr/bin/python3
# Gets distance and conservation information for motifs
# from other organisms.

import subprocess

s = "../../../../python/motif/distOnly.py"

fimoBase = "/media/jburdick/disk2/jburdick/fimo/cb3_"

outputDir = "/media/jburdick/disk2/jburdick/distAndConservation_cb3/"

upstreamRegions = "../../../../data/seq/cb3/cb3_Ce_gene_upstream_3kb.bed"

subprocess.call([s, fimoBase + "Ce_1.02", outputDir + "Ce_1.02", upstreamRegions])
subprocess.call([s, fimoBase + "Dm_1.02", outputDir + "Dm_1.02", upstreamRegions])
subprocess.call([s, fimoBase + "Hs_1.02", outputDir + "Hs_1.02", upstreamRegions])
subprocess.call([s, fimoBase + "Mm_1.02", outputDir + "Mm_1.02", upstreamRegions])

