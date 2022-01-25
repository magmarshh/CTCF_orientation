# CTCF_orientation

CTCF_orientation is a python tool that produces a pie chart visualizing the CTCF orientation of chromatin loops in a given loop file. 

## Requirements

- python 3.7
- pandas 
- matplotlib


## Usage

```{bash echo=FALSE}
ctcf_orientation.py -l <loopfile> -m <motiffile> -o <outputfile>
```
### Parameters

- -l; loop file: chromatin loop file in bedpe format (see [bedtools documentation](https://bedtools.readthedocs.io/en/latest/content/general-usage.html) for more information), with or without a header row. 
- -m; motif file: CTCF motif file in the format "chromosome start end name strength orientation pvalue qvalue sequence", although only "chromosome start end" will be used, with the other columns optional. 
- -o; output file: Desired file name for the outputted piechart with extensions available in matplotlib.
- h; help: Helpful information for usage. 

## Output

CTCF_orientation produces a pie chart showing the distribution of CTCF motifs in the chromatin loop dataset provided. CTCF motifs will be binned into 5 categories: convergent, divergent, tandem, single, and none depending on which CTCF motifs the first anchor of the loop intersects with and which the second does. For example, if the first anchor contains a CTCF site but the second doesn't, that loop is designated as containing a 'single' CTCF site. If both anchors contain CTCF sites on the same strand, the loop contains a tandem CTCF motif, etc. 

## Example 

```{bash echo=FALSE}
ctcf_orientation.py -l ./test_data/data/grubert_hg19_gm12878_original_loops.bed -m ./test_data/CTCF.bed -o ./test_data/gm12878_hg19_grubert_loops_ctcf_piechart.svg
```

## Thanks 
- [CTCF](https://github.com/mdozmorov/CTCF), Source for CTCF motif files and ideas. 






