# CTCF_orientation

CTCF_orientation is a python tool that produces a pie chart visualizing the CTCF orientation distribution of chromatin loops in a given loop file. 

## Requirements

- python 3.7
- pandas 
- matplotlib


## Usage

```{bash echo=FALSE}
ctcf_orientation.py -l <loopfile> -m <motiffile> -o <outputfile>
```
### Parameters

- -l; loop file: Pairedbed file of loops with the format "chromosome1 anchor1_start anchor1_end chromosome2 anchor2_start anchor2_end"
- -m; motif file: CTCF motif file in the format "chromosome start end name strength orientation pvalue qvalue sequence", although only "chromosome start end" will be used, with the other columns optional. 
- - -o; output file: Desired file name for the outputted piechart with extensions available in matplotlib.
- h; help: Helpful information for usage. 

## Example 

```{bash echo=FALSE}
ctcf_orientation.py -l ./test_data/data/grubert_hg19_gm12878_original_loops.bed -m ./test_data/CTCF.bed -o ./test_data/gm12878_hg19_grubert_loops_ctcf_piechart.svg
```

## Thanks 
- [CTCF](https://github.com/mdozmorov/CTCF), Source for CTCF motif files and ideas. 
- [CTCF_orientation](https://github.com/magmarshh/CTCF_orientation), source of core program and root of all code






