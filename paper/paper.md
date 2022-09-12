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

`Wakeflow` is a Python package for generating models of the perturbations 
induced by planets embedded in gaseous circumstellar disks. These perturbations
take the form of a spiral shock wave [@Ogilvie:2002], and are often called 
a "planet wake" in analogy with that produced by a boat in a lake.

Detecting newly formed planets embedded in their disk is a challenging problem 
in the field of planet formation. A major area of progress in recent years is
the detection of planets by the gravitationally induced disturbance in their
host disks. This disturbance, caused by the planet wake, manifests as a deviation 
in velocity from the bulk flow which may be measured through the Doppler shift of 
molecular lines [e.g. @Perez:2015; @Pinte:2018]. Such kinematic observations have 
been accurately reproduced through 3D fluid simulations of the planet-disk interaction, 
allowing for the inference of planet and disk properties [@Pinte:2018; @Pinte:2019]. 
However, such studies are computationally very expensive.

`Wakeflow` eases this computational cost by applying the theory of planet wake 
generation and propagation [@Goldreich:1979; @Goodman:2001; @Rafikov:2002; @Bollati:2021]
to create semi-analytic models of planet wakes. `Wakeflow` models are readily created in 
less than a minute on a modern laptop, as opposed to the hours of supercomputer time needed 
for 3D hydrodynamical simulations. The relatively low computational cost of `Wakeflow`
means that researchers can get an idea of whether planet-disk interactions
can explain their observations, and the disk and planet parameters needed, before
spending computer time on more detailed simulations. In addition, `Wakeflow`
allows researchers to create these models without needing detailed theoretical 
knowledge. 

Synthetic observations of disk models are common in the literature, produced
through radiation transfer codes such as `MCFOST` [@Pinte:2006; @Pinte:2009]. `Wakeflow` 
is capable of interfacing with `MCFOST` in order to create synthetic observations of the
semi-analytic models for direct comparison with real continuum or line emission
observations.

`Wakeflow` is partially adapted from a previous Python code also written by us called
`Analytical_Kinks` [@Analytical_Kinks]. `Wakeflow` is intended to be a more complete, versatile 
and easy to use version of that code, and it obeys standard Python packaging conventions.
In addition, `Wakeflow` can directly interface with `MCFOST` while `Analytical_Kinks` cannot.
At the time of writing, no other open source software packages exist to generate the perturbations
induced by an embedded planet in a circumstellar disk using the semi-analytic theory
of planet wakes.

Existing scientific publications focusing on detecting the kinematic signatures of
planets that have used `Wakeflow` or its predecessor `Analytical_Kinks` include
@Bollati:2021, @Calcino:2022 and @Teague:2022.

# Acknowledgements

`Wakeflow` relies on the following scientific Python packages: `NumPy` [@Harris:2020], 
`Matplotlib` [@Hunter:2007], `SciPy` [@Virtanen:2020] and `Astropy` [@Astropy:2022].

# References