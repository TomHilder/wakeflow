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
    orcid: 0000-0003-4679-4072
    affiliation: 2
  - name: Francesco Bollati
    affiliation: 3
  - name: Jacob Vandenberg
    orcid: 0000-0002-4055-9728
    affiliation: 4
affiliations:
 - name: School of Physics and Astronomy, Monash University, Australia
   index: 1
 - name: Dipartimento di Fisica, Università degli Studi di Milano, Italy
   index: 2
 - name: Dipartimento di Scienza e Alta Tecnologia, Università degli Studi dell'Insubria, Italy
   index: 3
 - name: School of Mathematics, Monash University, Australia
   index: 4
date: 19 September 2022
bibliography: paper.bib
---

# Summary

`Wakeflow` is a Python package for generating semi-analytic models of the perturbations 
induced by planets embedded in gaseous circumstellar disks. These perturbations
take the form of a spiral shock wave [@Ogilvie:2002], and are often called 
a "planet wake" in analogy with that produced by a boat in a lake.
Using `Wakeflow`, users may calculate the perturbed density and velocity fields of the gas in the disk.
These may be used with radiation transfer codes to generate synthetic observations
that map both the gas distribution and the gas kinematics. 
Comparison with real observations, such as from molecular line emission taken with the 
Attacama Large Millimetre Array, allows researchers to infer the properties of potential planets
as well as the disk itself.

# Statement of need

Detecting newly formed planets embedded in their disk is a challenging problem 
in the field of planet formation. A major area of progress in recent years is
the detection of planets by the gravitationally induced disturbance in their
host disks. This disturbance, caused by the planet wake, manifests as a deviation 
in velocity from the bulk flow which may be measured through the Doppler shift of 
molecular lines [e.g. @Perez:2015; @Pinte:2018]. Such kinematic observations have 
been accurately reproduced through 3D fluid simulations of the planet-disk interaction, 
allowing for the inference of planet and disk properties [@Pinte:2018; @Pinte:2019]. 
However, these studies are computationally expensive.

`Wakeflow` eases this computational cost by applying the theory of planet wake 
generation and propagation [@Goldreich:1979; @Goodman:2001; @Rafikov:2002; @Bollati:2021]
to create semi-analytic models of planet wakes. `Wakeflow` models are readily created in 
less than a minute on a modern laptop, as opposed to the hours of supercomputer time needed 
for 3D hydrodynamical simulations. The relatively low computational cost of `Wakeflow`
means that researchers can get an idea of whether planet-disk interactions
can explain their observations, and the disk and planet parameters needed, before
spending computer time on more detailed simulations.

`Wakeflow` can interface with the radiative transfer code `MCFOST` [@Pinte:2006; @Pinte:2009] in order to create synthetic observations of the
semi-analytic models for direct comparison with observed continuum or line emission.

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
