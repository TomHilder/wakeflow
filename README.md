<div id="top"></div>

<!-- PROJECT SHIELDS -->
[![PyPI Package][pypi-shield]][pypi-url]
[![Python][python-shield]][python-url]
[![License][license-shield]][license-url]
[![JOSS][JOSS-shield]][JOSS-url]
[![Docs][docs-status-shield]][docs-status-url]



<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/TomHilder/wakeflow">
    <img src="https://github.com/TomHilder/wakeflow/blob/main/logo.png?raw=true" alt="Logo" width="500" height="500">
  </a>



<!--  <h3 align="center">Wakeflow</h3> -->

  <p align="center">
    Generate and manipulate semi-analytic models of planet wakes
    <br />
    <br />
    <a href="https://wakeflow.readthedocs.io/en/latest/tutorials/quickstart.html"><strong>Quickstart tutorial »</strong></a>
    <br />
    <a href="https://wakeflow.readthedocs.io/en/latest/index.html"><strong>Documentation »</strong></a>
    <br />
  </p>
</div>


<!-- ABOUT THE PROJECT -->
# Overview

`wakeflow` is a Python package primarily for calculating tidally-induced perturbations resulting from a planet embedded in a gas disk. It is an implementation of both the linear theory for planet wake generation ([Goldreich and Tremaine 1979](https://ui.adsabs.harvard.edu/abs/1979ApJ...233..857G)) and the non-linear theory of wake propagation ([Rafikov 2002](https://ui.adsabs.harvard.edu/abs/2002ApJ...569..997R/abstract)) in 2D. `wakeflow` lets you generate these models by specifying disk and system properties as typically parameterised in the planet formation literature. It also contains additional tools allowing you to:
* Visualise your results
* Create 3D models under some assumptions
* Interface directly with the radiative transfer code [MCFOST](https://github.com/cpinte/mcfost) to generate synthetic images of these models
* (Planned) Rotate and project your models to create line-of-sight maps of velocity perturbations at some emitting layer
* (Planned) Create analytic predictions for peak velocity maps as found in [Calcino et al. 2022](https://ui.adsabs.harvard.edu/abs/2022ApJ...929L..25C/abstract)

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Installation

`wakeflow` may be most easily installed from the Python Package Index ([PyPI](https://pypi.org/project/wakeflow/)), or can also be installed from the [GitHub repository](https://github.com/TomHilder/wakeflow) if you wish make contributions. Dependencies for `wakeflow` consist mostly of standard python libraries. We recommend using a package manager such as [Anaconda](https://www.anaconda.com/products/distribution) to make your life easier, but this is not required.



### Dependencies

Python packages:

* `numpy`
* `matplotlib`
* `astropy`
* `scipy`
* `setuptools`
* `pyyaml`
* `tqdm`
* `pytest`
* [`pymcfost`](https://github.com/cpinte/pymcfost) (if interfacing with [MCFOST](https://github.com/cpinte/mcfost))



### PyPI (pip)

The easiest way to install `wakeflow` is via [PyPI](https://pypi.org/project/wakeflow/), using `pip`:
```sh
pip install wakeflow
```
that's it!



### From source (GitHub)

If you want to contribute to, or modify `wakeflow`, you should install it from the [GitHub repository](https://github.com/TomHilder/wakeflow). Simply fork the repo using the button in the top right, and then clone it:
```sh
git clone https://github.com/<replace-by-your-username>/wakeflow.git
```
Alternatively, you may install from source even if you do not with to edit `wakeflow`, in which case I would recommend skipping the fork and simply cloning the repo directly:
```sh
git clone https://github.com/TomHilder/wakeflow.git
```
Navigate to the directory it is installed in:
```sh
cd wakeflow
```
You can verify that you are in the correct directory by checking that you see `README.md` when you run:
```sh
ls
```
Now we use `pip` to create a local and _editable_ install of `wakeflow`:
```sh
python -m pip install -e .
```
Do not forget the dot (.) in the above command, as it tells `pip` to look in the current working directory (where `wakeflow` is). The advantage of installing this way is that it places a _link_ to `wakeflow` in your `site-packages` folder instead of _moving_ it there. Now when you edit the code in `wakeflow/src/wakeflow/` it will edit your installation!

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

Please refer to the [Quickstart tutorial](https://wakeflow.readthedocs.io/en/latest/tutorials/quickstart.html) for the most typical usage of `wakeflow` including generating models and reading the results. Additional examples of more advanced usage can be found in the [Documentation](https://wakeflow.readthedocs.io/en/latest/index.html).

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- CONTRIBUTING -->
## Contributing

Contributions to `wakeflow` are welcome. If you would like to implement a new feature, please:

1. Install using the above installation from source instructions
2. Create your Feature Branch (`git checkout -b feature/NewFeature`)
3. Commit your Changes (`git commit -m 'Added some NewFeature'`)
4. Push to the Branch (`git push origin feature/NewFeature`)
5. Open a Pull Request

If you have a suggestion that would improve `wakeflow` but do not have the time or means to implement it yourself, please simply open an issue with the tag "enhancement". If you would like to report a bug, please open an issue with the tag "bug".

Don't forget to give the project a star!

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- CITATION -->
## Citing

Please cite [Hilder et al. (in prep)](https://example.com) in any work where `wakeflow` has been used. Please contact us if `wakeflow` is useful to you, we welcome any collaboration opportunities.

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

If you are having difficulties installing or using `wakeflow`, please let us know! We are happy to answer any questions or provide assistance. 

Thomas Hilder - thil0004@student.monash.edu

Project Link: [https://github.com/TomHilder/wakeflow](https://github.com/TomHilder/wakeflow)

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

`wakeflow` is based on the semi-analytic theory of planets wakes described in [Rafikov (2002)](https://iopscience.iop.org/article/10.1086/339399) and [Bollati et al. (2021)](https://academic.oup.com/mnras/article-abstract/504/4/5444/6255419). The code is partially adapted from `analytical kinks` which was written by Francesco Bollati, Daniele Fasano and Thomas Hilder, and can be found [here](https://github.com/DanieleFasano/Analytical_Kinks).

Additional acknowledgements:
* [Img Shields](https://shields.io)
* [Best-README-Template](https://github.com/othneildrew/Best-README-Template)

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[pypi-shield]: https://img.shields.io/badge/PyPI-1.1.1-green?style=flat-square
[pypi-url]: https://pypi.org/project/wakeflow/
[python-shield]: https://img.shields.io/badge/Python-3.6%2B-orange?style=flat-square
[python-url]: https://www.python.org
[license-shield]: https://img.shields.io/badge/License-MIT-blue?style=flat-square
[license-url]: https://opensource.org/licenses/MIT
[JOSS-shield]: https://img.shields.io/badge/JOSS-coming%20soon-blueviolet?style=flat-square
[JOSS-url]: https://example.com
[docs-status-shield]: https://readthedocs.org/projects/wakeflow/badge/?version=latest&style=flat-square
[docs-status-url]: https://wakeflow.readthedocs.io/en/latest/?badge=latest