import os
from pathlib import Path

import pandas as pd


def get_chromsizes(species, chromsizes=None):
    """Get chromsizes file for species."""
    # extract species name from available files
    csizes = Path("data/chromsizes") 
    snames = [sp.split(".")[0] for sp in os.listdir(csizes)]
    if chromsizes is not None:
        csdf = pd.read_csv(chromsizes, sep='\t', header=None)
        csdict = pd.Series(csdf[1].values, index=csdf[0]).to_dict() 
    elif species not in snames:
        raise Exception("Species not yet supported, please provide a chromsizes file.")
    else: 
        csfile = csizes / f"{species}.chrom.sizes"
        csdf = pd.read_csv(csfile, sep='\t', header=None)
        csdict = pd.Series(csdf[1].values, index=csdf[0]).to_dict() 
    return(csdict)

#get_chromsizes("equCab2", chromsizes="/home/groups/hoolock/u1/genomes/horse/EquCab2/equCab2.chrom.sizes")
