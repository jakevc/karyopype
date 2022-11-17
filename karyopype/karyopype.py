import pathlib
from typing import Optional
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.collections import BrokenBarHCollection
from matplotlib.lines import Line2D
from pkg_resources import resource_filename, resource_listdir


def filter_cannonical(df):
    """Retains only cannonical chromosomes in bed-like df."""
    df = df[df.loc[:, 0].str.match("chr[0-9|X|Y]+[a|b]?$")]
    return df


def list_species():
    data_files = resource_listdir(__name__, 'data/chromsizes/')
    splist = [sp.split('.')[0] for sp in data_files]
    print(splist)
    return(splist)


def get_chromsizes(species, chromsizes=None, cannonical=True):
    """Return chromsizes dict for speceis."""
    # check that species is a string
    if not isinstance(species, str):
        raise TypeError("Species name should be a string, e.g 'hg38'")

    # fetch the pkg chromsizes data
    data_files = resource_listdir(__name__, 'data/chromsizes/')
    snames = [sp.split('.')[0] for sp in data_files]

    if chromsizes is not None:
        csdf = pd.read_csv(chromsizes, sep='\t', header=None)
        # filter cannonical by default for now
        csdf = filter_cannonical(csdf)
        csdict = pd.Series(csdf[1].values, index=csdf[0]).to_dict()
    elif species not in snames:
        raise Exception(
            "Species not yet supported, please provide a chromsizes file.")
    else:
        csname = [s for s in data_files if species in s][0]
        csfile = resource_filename(__name__, f'data/chromsizes/{csname}')
        csdf = pd.read_csv(csfile, sep='\t', header=None)
        csdf = filter_cannonical(csdf)
        csdict = pd.Series(csdf[1].values, index=csdf[0]).to_dict()

    return(csdict)


def parse_regions(regions=None,sep:Optional='\t'):
    """Return regions file if any, skip is True if None."""
    # determine if there is a dataframe of other bed regions
    if (regions is not None
            and isinstance(regions, str)
            or isinstance(regions, pathlib.PosixPath)):
        regions = pd.read_csv(regions, sep=sep, header=None).iloc[:, 0:3]
        regions.columns = ["chrom", "start", "end"]
        skip = False
    elif isinstance(regions, pd.DataFrame):
        regions = regions.iloc[:, 0:3]
        # add expected columns to reigons
        regions.columns = ["chrom", "start", "end"]
        skip = False
    else:
        skip = True
        print("No additional regions, only showing chromsizes")
    return regions, skip


def add_chromsize_colors(sizes_dict, color):
    """Add matplotlib color to chromsizes dict."""
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


def add_regions(ax, chromsizes, regions=None,sep:Optional='\t'):
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
    
    # legend
    leg = []
    leg_lab = []

    if regions is None:
        print("No additional regions.")
    elif isinstance(regions, str):
        # add single region as green regions
        rdf, skip = parse_regions(regions,sep=sep)
        color = 'C1'
        if skip is False:
            print("Writing additional regions to axis.")
            rdf['colors'] = color

            # add legend parts
            leg.append(Line2D([0], [0], color=color, lw=4))
            leg_lab.append("regions1")

            for collection in chromosome_collections(
                    rdf, chrom_ybase, chrom_height,
                    edgecolor=color):
                ax.add_collection(collection)

    else:
        for i, r in enumerate(regions):
            color = f'C{i}'

            # add legend element per regions
            leg_lab.append(f"regions{i+1}")
            leg.append(Line2D([0], [0], color=color, lw=4))

            # determine if there is a dataframe of other bed regions
            rdf, skip = parse_regions(r)
            if skip is False:
                rdf['colors'] = color
                for collection in chromosome_collections(
                        rdf, chrom_ybase, chrom_height,
                        edgecolor=color):
                    ax.add_collection(collection)

    # add to ax
    ax.set_yticks([chrom_centers[i] for i in chromosome_list])
    ax.set_yticklabels(chromosome_list)

    ax.legend(leg, leg_lab, loc=4)
    ax.axis("tight")
    return(ax)


def plot_karyopype(species, regions=None,
                   chromsizes=None, savefig=False,
                   figsize=(10, 7)):
    """
    Plot karyopype of the genome of interest,
    along with an extra set of genomic regions or a list of genomic regions.
        - species:
            The name of the species chromosomes to plot, e.g. 'hg38', 'nomLeu3'
        - chromsizes:
            Return regions file if any, skip is True if None.
            If not available in: karyopype.list_species().
        - regions:
            A dataframe, file or list with "chr, start, end" as first three columns.
        - savefig:
            Saves the chromosome plot to a file.
    """
    # call chromsizes
    chromsizes = get_chromsizes(species, chromsizes)

    # initialize a figure/
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
