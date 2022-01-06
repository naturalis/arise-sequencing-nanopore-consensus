# arise-sequencing-nanopore-consensus
## starting point
This repo should be considered as a follow-up of [arise-sequencing-nanopore](https://github.com/naturalis/arise-sequencing-nanopore) and assumes
data has been demultiplexed, either with *guppy* (prefered) or in *bash* (e.g. [demultiplex amplicons](https://github.com/naturalis/arise-sequencing-nanopore#demultiplex-amplicons-specimens-within-datasets-bash)).

## create 'nanopore' environment (Conda)
Requirement: conda (anaconda/miniconda)

`conda create -n nanopore python=3.6 pip`\
`conda activate nanopore`\
`conda install --yes -c conda-forge -c bioconda medaka==0.11.5 openblas==0.3.3 spoa racon minimap2`\
`pip install NGSpeciesID`\
`conda install -c bioconda prinseq cutadapt seqtk`

## consensus calling
The following command filters reads 25 nt above and below the median with Prinseq.\
Consensus calling is done with NGSpeciesID.\
Fasta headers and filenames were replaced using a [name_code_index file](https://github.com/naturalis/arise-sequencing-nanopore-consensus/blob/main/index_files/name_code_index_fungi.txt)\
Primers of renamed consensus sequences were trimmed using Cutadapt;\
sequence orientation was corrected (revcomp) if necessary.
From within the folder of demultiplexed fastqs run:

`./ont_consensus.sh $1 $2 $3`

$1 = name_code_index.txt, $2 = forward primer, $3 = revcomp of reverse primer\
[The ont_consensus.sh script is located here](https://github.com/naturalis/arise-sequencing-nanopore-consensus/tree/main/scripts/ont_consensus.sh)

The output file (out_trimmed.fasta) has been added as [out_trimmed_fungi.fasta](https://github.com/naturalis/arise-sequencing-nanopore-consensus/blob/main/metadata/out_trimmed_fungi.fasta) to the metadata folder.

## blast search consensus sequences
Blast searches were done against the fungi.fas reference dataset using blastn in Galaxy. All sequences had 100% coverage and the average identity was%. A [summary of the fungi blast results](https://github.com/naturalis/arise-sequencing-nanopore-consensus/blob/main/metadata/fungi_blast.md) suggests something might have gone wrong with the indices in the wetlab.

The output [fungi_blast_results.xlsx](https://docs.google.com/spreadsheets/d/1O1ae-nSbFfYfkCb8iElj5CIIvgTfNb5V/edit?usp=sharing&ouid=109237925768461347094&rtpof=true&sd=true) has been added to the Arise sequencing folder.
