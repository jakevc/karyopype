#!/bin/bash
# fetch chromsizes files frome UCSC, deposit in current directory.
for sp in hg38 mm10 nomLeu3 rheMac8 panTro4 gorGor3
    do wget http://hgdownload.cse.ucsc.edu/goldenPath/$sp/bigZips/$sp.chrom.sizes
done
