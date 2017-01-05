assimilr
=====

Two localized algorithms based on the EnKPF: the naive-LEnKPF and the block-LEnKPF.

## Description

Implementation of various localized algorithms based on the EnKPF (introduced in [Frei & Kunsch](http://biomet.oxfordjournals.org/content/100/4/781.short)). 
The new local algorithms are the naive-LEnKPF and block-LEnKPF. For details see [Robert & Kunsch](https://arxiv.org/abs/1605.05476), soon to appear in print. 

The algorithms are tested on the modified SWEQ model of [Wursch & Craig 2014](http://www.meteo.physik.uni-muenchen.de/dokuwiki_en/lib/exe/fetch.php?media=lscraig:herz:sw_model_submitted.pdf), 
here implemented in fortran90 for better performance. 

## Instructions

To install the package use devtools::install_github("robertsy/assimilr"). The folder example provides some basic usage of the main functions, 
once for a one-step assimilation in a one-dimensional Gaussian random field and then for more complicated cases with the modified
SWEQ model. 
