[![CircleCI](https://circleci.com/gh/jakevc/karyopype.svg?style=svg)](https://circleci.com/gh/jakevc/karyopype)

# karyopype

Karyopype is a simple chromosome plotting package in python allowing you to quickly visualize where a set of genomic regions, or multiple sets of genomic regions fall on the chromosomes.


```
import karyopype
kp = Karyopype("hg38")
kp.add_regions(regions_bed_df)
kp.plot_chromosomes()
```