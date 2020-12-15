# Annotathon

Methodology for the [Annotathon Project](http://annotathon.org/)

There is a description of the project on the site, available both in French and English. We encourage you to read it thouroughly before diving into this repository. Everything will naturally make sense.

Authors :

- MAGAÑA LOPEZ Gustavo

- BOUTARD Anthony

We want to maximise our workflow's reproducibility, first and foremost because this is a personal conviction, and secondly because the project explicitly requires us to. We have thus created this repo. 

Manually annotating metagenomic sequences is tiresome and error-prone when done completely by hand. Many different sites have to be consulted, and sequences have to be copied and pasted repeatedly. The graphical interfaces of these websites are full of parameters that control the [ORF](https://en.wikipedia.org/wiki/Open_reading_frame) finding procedure, to name the first step of the annotation workflow. This is problematic because one must take the time to carefully write down all parameters or take a screenshot for later logging. 

This is just the first problem. We believe that the process of annotating new sequences should focus on discussing the data, and looking into potential anomalies and/or discoveries. The bioinformatician's time and energy should prioritise analysis rather than book-keeping and mechanic repetition of interaction with websites.

__This document, and most of the project will be documented in a bilingually, both in French and in English (mainly). French because it is the language of our academic community. However, English is the _de facto lingua franca_ of science. Documenting the repository in English will broaden its international reach within the bioinformatics community.__



## R support

We decided to add a directory `R/` to leverage BioConductor's wonderfull tools. 

This includes a [customised version](https://github.com/gmagannaDevelop/msa-mod)  of the package `msa` to perform Multiple Sequence Alignment, using parallel processing to optimise gap parameters.



## Extra :

Take a look into the `misc` directory. There you will find various documents like a manually created changelog and a todo list.
Some caveats and checksums are also provided.
