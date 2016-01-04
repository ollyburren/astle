#  ASTLE

This repo contains scripts that have been developed to analyse Will Astle's GWAS for blood traits and combine with functional and promoter capture Hi-C data.
It requires CHIGP repo to be installed and access to correctly formatted support datasets is mandatory. This means that it is only of use to Olly Burren at the moment, however when the analysis is published I will revisit and make the data available so that others can rerun the analysis and tweak.

## Recipe

I assume that the data from tar ball has been extracted and sits in DATA dir in this repo.

Convert summary statistics into pvalues for inclusion.
``perl ./perl/computePPi.pl``
Combine CHiC, and other functional information to compute a gene score.
``perl ./perl/geneScores.pl``

## TODO

Make data available on publication
