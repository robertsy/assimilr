assimilr
=====

Two localized algorithms based on the EnKPF: the naive-LEnKPF and the block-LEnKPF.

## Description

Implementation of two new local algorithms based on the EnKPF  ([Frei & Kunsch 2013](http://biomet.oxfordjournals.org/content/100/4/781.short)):
the naive-LEnKPF and the block-LEnKPF. For details see [Robert & Kunsch 2017](http://www.tandfonline.com/doi/full/10.1080/16000870.2017.1282016). 

The algorithms are tested on the modified SWEQ model of [Wursch & Craig 2014](http://www.meteo.physik.uni-muenchen.de/dokuwiki_en/lib/exe/fetch.php?media=lscraig:herz:sw_model_submitted.pdf), 
here implemented in fortran90 for better performance. 

## Instructions

To install the package use devtools::install_github("robertsy/assimilr"). The folder example provides some basic usage of the main functions, 
once for a one-step assimilation in a one-dimensional Gaussian random field and then for more complicated cases with the modified
SWEQ model. 
