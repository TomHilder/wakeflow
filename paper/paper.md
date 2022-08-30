---
title: 'Wakeflow: A Python package for semi-analytic models of planetary wakes'
tags:
  - Python
  - astronomy
  - dynamics
  - protoplanetary disks
  - planet formation
authors:
  - name: Thomas Hilder
    orcid: 0000-0001-7641-5235
    corresponding: true
    affiliation: 1
  - name: Daniele Fasano
    orcid: 0000-0000-0000-0000
    affiliation: 2
  - name: Francesco Bollati
    orcid: 0000-0000-0000-0000
    affiliation: 3
affiliations:
 - name: School of Physics and Astronomy, Monash University, Australia
   index: 1
 - name: Dipartimento di Fisica, Università degli Studi di Milano, Italy
   index: 2
 - name: Dipartimento di Scienza e Alta Tecnologia, Università degli Studi dell'Insubria, Italy
   index: 3
date: 26 August 2022
bibliography: paper.bib
---

# Summary

`wakeflow` is a Python package for generating models of the perturbations 
induced by planets embedded in their natal gas disks. These perturbations
take the form of a spiral shock wave (cite Ogilvie), and are often called 
a "planet wake" in analogy with that produced by a boat in a lake.

Detecting newly formed planets embedded in their disk is a challenging problem 
in the field of planet formation. A major area of progress in recent years is
the detection of planets by the gravitationally induced disturbance in their
host disks. This disturbance manifests as a deviation in velocity from the bulk
flow which may be measured through the Doppler shift of molecular lines (citations). 
Such kinematic observations have been accurately reproduced through 3D fluid simulations
of the planet-disk interaction, allowing for the inference of planet and disk
properties (citations). However, such studies are computationally very expensive.

`wakeflow` eases this computational cost by applying the theory of planet wake 
generation and propagation (goodman, rafikov, bollati) to create semi-analytic 
models of planet wakes. `wakeflow` models are readily created in less than a 
minute on a modern laptop, as opposed to the often hours of supercomputer time needed 
for 3D hydrodynamical simulations. The relatively low computational cost of `wakeflow`
means that researchers can get an idea of whether planet-disk interactions
can explain their observations, and the disk and planet parameters needed, before
spending computer time on more detailed simulations. In addition, `wakeflow`
allows researchers to create these models without needing detailed theoretical 
knowledge. 

Synthetic observations of disk models are common in the literature, produced
through radiation transfer codes such as `MCFOST` (pinte). `wakeflow` is capable of 
interfacing with `MCFOST` in order to create synthetic observations of the
semi-analytic models for direct comparison with real continuum or line emission
observations.

# Statement of need



# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Acknowledgements



# References