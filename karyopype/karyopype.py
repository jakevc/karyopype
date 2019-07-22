import os
from pathlib import Path
import pathlib

# import matplotlib._color_data as mcd
# import matplotlib.patches as mpatch
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.collections import BrokenBarHCollection


def filter_cannonical(df):
    """Retains only cannonical chromosomes in first column."""
    df = df[df.loc[:, 0].str.match("chr[0-9|X|Y]+[a|b]?$")]
    return df


def get_chromsizes(species, chromsizes=None, cannonical=True):
    """Get chromsizes file for species."""
    # check that species is a string
    if not isinstance(species, str):
        raise TypeError("Species name should be a string, e.g 'hg38'")

    # extract species name from available files
    csizes = Path("data/chromsizes/")
    snames = [sp.split(".")[0] for sp in os.listdir(csizes)]
    if chromsizes is not None:
        csdf = pd.read_csv(chromsizes, sep='\t', header=None)
        # filter cannonical by default for now
        csdf = filter_cannonical(csdf)
        csdict = pd.Series(csdf[1].values, index=csdf[0]).to_dict()
    elif species not in snames:
        raise Exception("Species not yet supported,\
            please provide a chromsizes file.")
    else:
        csfile = csizes / f"{species}.chrom.sizes"
        csdf = pd.read_csv(csfile, sep='\t', header=None)
        csdf = filter_cannonical(csdf)
        csdict = pd.Series(csdf[1].values, index=csdf[0]).to_dict()

    return(csdict)


def parse_regions(regions=None):
    # determine if there is a dataframe of other bed regions
    if (regions is not None
            and isinstance(regions, str)
            or isinstance(regions, pathlib.PosixPath)):
        regions = pd.read_csv(regions, sep=' ', header=None)
        regions.columns = ["chrom", "start", "end"]
        skip = False
    elif isinstance(regions, pd.DataFrame):
        regions = regions
        # add expected columns to reigons
        regions.columns = ["chrom", "start", "end"]
        skip = False
    else:
        skip = True
        print("No additional regions, only showing chromsizes")
    return regions, skip


def add_chromsize_colors(sizes_dict, color):
    """Add matplotlib color to chromsizes dict"""
    for k in sizes_dict:
        sizes_dict[k] = (sizes_dict[k], color)
    return(sizes_dict)


# the following functions are heavily inspired by the follwoing gist:
# https://gist.github.com/daler/c98fc410282d7570efc3
def chromosome_collections(df, y_positions, height, **kwargs):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    for chrom, group in df.groupby('chrom'):
        yrange = (y_positions[chrom], height)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=group['colors'], **kwargs)
    if del_width:
        del df['width']


def add_regions(ax, chromsizes, regions=None):
    """Add Chromosome breakpoints plot to axis."""

    chromsize_colors = add_chromsize_colors(chromsizes, color="whitesmoke")
    chrom_height = 1  # Height of each ideogram
    chrom_spacing = 1  # Spacing between consecutive ideograms

    gene_height = 0.4  # Height of the gene track < `chrom_spacing`
    gene_padding = 0.1  # Padding between the top of a gene track ideogram

    # Decide which chromosomes to use
    chromosome_list = list(chromsizes.keys())

    # Keep track of the y positions genes for each chromosome
    # and center of each ideogram (which is where we'll put the ytick labels)
    ybase = 0
    chrom_ybase = {}
    gene_ybase = {}
    chrom_centers = {}

    # Items in the beginning of `chromosome_list` will
    # appear at the top of the plot
    for chrom in chromosome_list[::-1]:
        chrom_ybase[chrom] = ybase
        chrom_centers[chrom] = ybase + chrom_height / 2.
        gene_ybase[chrom] = ybase - gene_height - gene_padding
        ybase += chrom_height + chrom_spacing

    # add chromsizes to ax
    ideo = pd.DataFrame(chromsize_colors).transpose().reset_index()
    ideo['start'] = 0
    ideo.columns = ["chrom", "end", "colors", "start"]

    # Add a new column for width
    ideo['width'] = ideo.end - ideo.start

    # Add a collection of the query chromsizes
    for collection in chromosome_collections(
            ideo, chrom_ybase, chrom_height, edgecolor='k'):
        ax.add_collection(collection)

    # determine if there is a dataframe of other bed regions
    rdf, skip = parse_regions(regions)

    if not skip:
        print("Writing additional regions to axis.")
        rdf['colors'] = 'red'
        for collection in chromosome_collections(
                rdf, chrom_ybase, chrom_height, edgecolor='red'):
            ax.add_collection(collection)

    # add to ax
    ax.set_yticks([chrom_centers[i] for i in chromosome_list])
    ax.set_yticklabels(chromosome_list)
    ax.axis("tight")
    return(ax)


def plot_karyopype(species, regions=None, savefig=False):
    """
    Plot karyopype of the genome of interest,
    along with an extra set of genomic regions or a list of genomic regions.
        - species:
            The name of the species chromosomes to plot, e.g. 'hg38', 'nomLeu3'
        - regions:
            A bedlike file or dataframe with at least chr, start, end.
        - savefig:
            Saves the chromosome plot to a file.
    """
    # call chromsizes
    chromsizes = get_chromsizes(species)

    # initialize a figure/
    figsize = (7, 5)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # add regions to the axis
    add_regions(ax=ax, chromsizes=chromsizes, regions=regions)

    # xtick formatting: position in Mbps
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, n: int(x/1e6)))

    # labels
    ax.set_title(f"{species} Karyopype", fontsize=14)
    plt.xlabel('Chromosome Position (Mbp)', fontsize=14)
    plt.ylabel(f'{species} Chromosome', fontsize=14)
    if savefig is True:
        plt.savefig(f'{species}_karyopype.png')

    return(plt)
