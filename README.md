# split-oligo-designer

## Overview

This repository contains the code to design the oligonucleotides for mFISH3D.

The code is heavily relying on [OligoMiner](http://dx.doi.org/10.1073/pnas.1714530115) tool. 



## Software prerequisites
Tested on Ubuntu 20.04 LTS with the following version of software.

- [BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (2.13.0)
- Python 3.8

Make sure to set the path to the BLAST.

## Installing dependencies

1. Make sure you have [conda](https://docs.conda.io/en/latest/miniconda.html) installed.
2. Create a new environement (python=3.8), and install the dependencies.

```
$ git clone https://github.com/tatz-murakami/split-oligo-designer.git
$ cd split-oligo-designer
$ conda env create -f environment.yml
$ conda activate oligo
```

## Database and template sequence

### Make database
1. Download and make the database for blast. The database for the cDNA  can be found at [Ensembl](http://www.ensembl.org/info/data/ftp/index.html.)
2. Generate the blast database according to the instruction of [NCBI](https://www.ncbi.nlm.nih.gov/books/NBK569841/). 
	example: `./data/human_transcriptome`
	
### Download fasta template
Download the fasta template from [NCBI](https://www.ncbi.nlm.nih.gov/) or else.
	example: `./data/slc17a7/slc17a7_hs.fasta`

## Example usage
The jupyter notebook is included `OligoDesign.ipynb`.

## Parameters
The example of the parameters are shown below.
```python
mFISH3D_param = {
    'fasta':'/home/tmurakami/src/split-oligo-designer/data/slc17a7/slc17a7_hs.fasta',
    'database':'/home/tmurakami/src/split-oligo-designer/data/human_transcriptome/human_transcriptome_db',
    'minimum_offtarget_gap':100,
    'hcr_seqs':{
            'seq_even_l':'GAGGAGGGCAGCAAACGGaa',
            'seq_odd_r':'atGAAGAGTCTTCCTTTACG',
            'seq_even_r':'',
            'seq_odd_l':''
        },
    'self_remove': True # set this True if your template sequence appears in database.
}
```

`minimum_offtarget_gap`
If the gap between two non-specific binding is more than minimum_offtarget_gap, the pair is not regarded to cause a off-target signal. Recommended value: 100.

`hcr_seqs`
The sequences of HCR fragments. The example design can be found in `./oligodesigner/parameters.py`

`self_remove` 
Set self_remove True if you want to remove the template sequence from
off-target analysis. Otherwise, the gene of your interest could be regarded as an off-target product.
Turn this to False when your template is not found in database (e.g. GFP).



The code requires the parameters for OligoMiner. The example is below.
```python
oligominer_param = {
    'l':20,
    'L':20,
    'gcPercent':25,
    'GCPercent':75,
    'tm':20,
    'TM':100,
    'X':'AAAAAA,TTTTTT,CCCCCC,GGGGGG',
    'sal':390,
    'form':30,
    'sp':1,
    'concA':25,
    'concB':25,
    'headerVal':None,
    'bedVal':False,
    'OverlapModeVal':False,
    'verbocity':False,
    'reportVal':True,
    'debugVal':False,
    'metaVal':False,
    'outNameVal':None,
    'nn_table':'DNA_NN3'
}
```


## Citation

> Murakami and Heintz. Multiplexed and scalable cellular phenotyping toward the standardized three-dimensional human neuroanatomy. bioRxiv (2022) https://doi.org/10.1101/2022.11.23.517711 

> Beliveau, B.J., Kishi, J.Y., Sasaki, H.M. et al. OligoMiner provides a rapid, flexible environment for the design of genome-scale oligonucleotide in situ hybridization probes. PNAS (2018) https://doi.org/10.1073/pnas.1714530115
	


## License

We provide this open source software without any warranty under the [MIT license](https://opensource.org/licenses/MIT).
