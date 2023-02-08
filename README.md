<div id="top"></div>

<!-- PROJECT SHIELDS -->
[![PyPI Package][pypi-shield]][pypi-url]
[![Python][python-shield]][python-url]
[![License][license-shield]][license-url]
[![JOSS][JOSS-shield]][JOSS-url]


![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/TomHilder/wakeflow/Tests.yml?label=tests&style=flat-square)
[![Docs][docs-status-shield]][docs-status-url]


<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/TomHilder/wakeflow">
    <img src="https://github.com/TomHilder/wakeflow/blob/master/logo.png?raw=true" alt="Logo" width="500" height="500">
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
* Interface directly with the radiative transfer code [`MCFOST`](https://github.com/cpinte/mcfost) to generate synthetic images of these models
* (Planned) Rotate and project your models to create line-of-sight maps of velocity perturbations at some emitting layer
* (Planned) Create analytic predictions for peak velocity maps as found in [Calcino et al. 2022](https://ui.adsabs.harvard.edu/abs/2022ApJ...929L..25C/abstract)


`wakeflow` is intended to allow both theorists and observers to easily generate models of the interaction between disks and embedded planets, instead of having to run expensive fluid simulations. In particular, `wakeflow` allows researchers to easily test whether a planet can explain kinematic perturbations observed in some set of disk observations, so-called _velocity kinks_ (see for example [Pinte et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...860L..13P/abstract)). `wakeflow` therefore also allows for a fast exploration of disk and planet parameters in order to determine those most likely to recreate observations, before resources are spent on 3D simulations. In addition, `wakeflow` models may be used with [`MCFOST`](https://github.com/cpinte/mcfost) to create synthetic images that may be compared directly with observations. 

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Installation

`wakeflow` may be most easily installed from the Python Package Index ([PyPI](https://pypi.org/project/wakeflow/)), or can also be installed from the [GitHub repository](https://github.com/TomHilder/wakeflow) if you wish to make contributions. Dependencies for `wakeflow` consist mostly of standard python libraries. We recommend using a package manager such as [Anaconda](https://www.anaconda.com/products/distribution) to make your life easier, but this is not required.




### PyPI (pip)

The easiest way to install `wakeflow` is via [PyPI](https://pypi.org/project/wakeflow/), using `pip`:
```sh
pip install wakeflow
```
that's it!



### From source (GitHub)

If you want to contribute to, or modify `wakeflow`, you should install it from the [GitHub repository](https://github.com/TomHilder/wakeflow). After installing the dependencies (see below), simply fork the repo using the button in the top right, and then clone it:
```sh
git clone https://github.com/<replace-by-your-username>/wakeflow.git
```
Alternatively, you may install from source even if you do not wish to edit `wakeflow`, in which case I would recommend skipping the fork and simply cloning the repo directly:
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



### Dependencies

Python packages:

* `numpy`
* `matplotlib`
* `astropy`
* `scipy`
* `setuptools`
* `pyyaml`
* `tqdm`
* `pytest` (optional)
* `pytest-cov` (optional)
* [`pymcfost`](https://github.com/cpinte/pymcfost) (optional, only if interfacing with [MCFOST](https://github.com/cpinte/mcfost))

If you install `wakeflow` using `pip` then the required dependencies will be automatically installed.



<!-- USAGE EXAMPLES -->
## Usage

Please refer to the [Quickstart tutorial](https://wakeflow.readthedocs.io/en/latest/tutorials/quickstart.html) for the most typical usage of `wakeflow` including generating models and reading the results. Additional examples of more advanced usage can be found in the [Documentation](https://wakeflow.readthedocs.io/en/latest/index.html).

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- USAGE EXAMPLES -->
## Testing

`wakeflow` is automatically unit-tested on Github using Actions and [`tox`](https://github.com/tox-dev/tox). If you have installed `wakeflow` from source, you may run a local test on your machine provided that you have `pytest` and `pytest-cov` installed. Simply navigate to your installation directory and run:
```sh
pytest
```


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
## Getting Help

If you are experiencing issues with `wakeflow`, please try the following:

1. Check the [documentation](https://wakeflow.readthedocs.io/en/latest/index.html) to see if it may be easily resolved
2. Open an [issue](https://github.com/TomHilder/wakeflow/issues) on the [repository](https://github.com/TomHilder/wakeflow)
3. Feel free to contact us directly using the details below


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
* [Tox](https://github.com/tox-dev/tox)
* [Tox GH Actions](https://github.com/ymyzk/tox-gh-actions)

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[pypi-shield]: https://img.shields.io/badge/pypi-1.2.0-green?style=flat-square
[pypi-url]: https://pypi.org/project/wakeflow/
[python-shield]: https://img.shields.io/badge/python-3.6%2B-orange?style=flat-square
[python-url]: https://www.python.org
[license-shield]: https://img.shields.io/badge/license-MIT-blue?style=flat-square
[license-url]: https://opensource.org/licenses/MIT
[JOSS-shield]: https://img.shields.io/badge/JOSS-under_review-yellow?style=flat-square
[JOSS-url]: https://joss.theoj.org/papers/22d35f9b0bb35a8df4dab85b0b6f4eb7
[test-status-shield]: https://img.shields.io/github/actions/workflow/status/TomHilder/wakeflow/Tests.yml?label=tests&style=flat-square
[test-status-url]: https://github.com/TomHilder/wakeflow/actions/workflows/Tests.yml
[docs-status-shield]: https://readthedocs.org/projects/wakeflow/badge/?version=latest&style=flat-square
[docs-status-url]: https://wakeflow.readthedocs.io/en/latest/?badge=latest
