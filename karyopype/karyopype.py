import os
from pathlib import Path

import matplotlib._color_data as mcd
import matplotlib.patches as mpatch
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.collections import BrokenBarHCollection

import genomes as gn


class Karyopype:
    """karyopype the given UCSC genome."""

    def __init__(self, species):
        
        self.species = species

        # check that species is a string
        if not isinstance(self.species, str):
            raise TypeError("Species name should be a string, e.g 'hg38'")

        self.chromsizes = gn.get_chromsizes(self.species)

    def add_chromsize_colors(self, sizes_dict, color):
        """Add matplotlib color to chromsizes dict"""
        for k in sizes_dict:
            sizes_dict[k]=(sizes_dict[k], color)
        return(sizes_dict)

    # the following functions are heavily inspired by the follwoing gist: 
    # https://gist.github.com/daler/c98fc410282d7570efc3
    def chromosome_collections(self, df, y_positions, height, **kwargs):
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

    def add_regions(self, ax, regions=None):
        """Add Chromosome breakpoints plot to axis."""

        chromsize_colors = self.add_chromsize_colors(self.chromsizes, color="whitesmoke")
        chrom_height = 1 # Height of each ideogram
        chrom_spacing = 1 # Spacing between consecutive ideograms
        
        gene_height = 0.4 # Height of the gene track < `chrom_spacing` 
        gene_padding = 0.1 # Padding between the top of a gene track ideogram
 
        # Width, height (in inches)
        figsize = (10, 8)

        # Decide which chromosomes to use
        chromosome_list = list(self.chromsizes.keys())

        # Keep track of the y positions for ideograms and genes for each chromosome,
        # and the center of each ideogram (which is where we'll put the ytick labels)
        ybase = 0
        chrom_ybase = {}
        gene_ybase = {}
        chrom_centers = {}

        # Iterate in reverse so that items in the beginning of `chromosome_list` will
        # appear at the top of the plot
        for chrom in chromosome_list[::-1]:
            chrom_ybase[chrom] = ybase
            chrom_centers[chrom] = ybase + chrom_height / 2.
            gene_ybase[chrom] = ybase - gene_height - gene_padding
            ybase += chrom_height + chrom_spacing

        ##### ADD CHROMSIZES #####

        # add chromsizes to ax 
        ideo = pd.DataFrame(chromsize_colors).transpose().reset_index()
        ideo['start'] = 0
        ideo.columns = ["chrom", "end", "colors", "start"]

        # Filter out chromosomes not in our list
        ideo = ideo[ideo.chrom.apply(lambda x: x in chromosome_list)]
        # Add a new column for width
        ideo['width'] = ideo.end - ideo.start

        # Add a collection of the query chromsizes
        for collection in self.chromosome_collections(ideo, chrom_ybase, chrom_height, edgecolor='k'):
            ax.add_collection(collection)

        # ### Add extra regions ####
        # for collection in self.chromosome_collections(self.bos, chrom_ybase, chrom_height, edgecolor='green'):
        #     ax.add_collection(collection)
            
        ### Add any regions
        
        # determine if there is a dataframe of other bed regions
        if isinstance(regions, str):
            regions = pd.read_csv(regions, sep='\t', header=None)
            regions.columns=["chrom","start","end","width"]
            skip = False
        elif isinstance(regions, pd.DataFrame):
            regions = regions
            # add expected columns to reigons
            regions.columns=["chrom","start","end","width"]
            skip = False
        else:
            skip=True
        
        if not skip:
            print("Writing additional regions to axis.")
            print(regions.columns)
            regions = regions.assign(colors = 'red')
            for collection in self.chromosome_collections(regions,
                                                          chrom_ybase,
                                                          chrom_height,
                                                          edgecolor='red'):
                ax.add_collection(collection)

        # add to ax
        ax.set_yticks([chrom_centers[i] for i in chromosome_list])
        ax.set_yticklabels(chromosome_list)
        ax.axis("tight")
        return(ax)

    def plot_chromosomes(self, regions=None, savefig=False):
        """
        Plot query-target breaks of synteny on query chromosomes
            - regions:
                A bedlike file or dataframe with at least chromosome, start, end.
            - savefig:
                Saves the chromosome plot to a file.
        """
    
        # initialize a figure/
        figsize = (10, 8)
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        
        # add regions to the axis
        self.add_regions(ax=ax, regions=None)
        
        # xtick formatting: position in Mbps
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x,n: int(x/1e6)))
        
        # labels
        ax.set_title(f"{self.species} Karyopype",fontsize=14)
        plt.xlabel('Chromosome Position (Mbp)', fontsize=14)
        plt.ylabel(f'{self.species} Chromosome', fontsize=14)
        plt.show()
        
        if savefig:
            plt.savefig(f"{self.species}karyopype_.png")
        

kp = Karyopype("hg19")
print(kp.chromsizes)
kp.plot_chromosomes()
