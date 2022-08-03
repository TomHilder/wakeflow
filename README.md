<div id="top"></div>

<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links

[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]
-->

[![PyPI Package][pypi-shield]][pypi-url]
[![Python][python-shield]][python-url]
[![License][license-shield]][license-url]
[![JOSS][JOSS-shield]][JOSS-url]

<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/TomHilder/wakeflow">
    <img src="https://github.com/TomHilder/wakeflow/blob/main/logo.png?raw=true" alt="Logo" width="500" height="500">
  </a>

<!--  <h3 align="center">Wakeflow</h3> -->

  <p align="center">
    Generate and manipulate semi-analytic models for planet wakes
    <br />
    <br />
    <a href="https://github.com/TomHilder/wakeflow"><strong>Quickstart tutorial »</strong></a>
    <br />
    <a href="https://github.com/TomHilder/wakeflow"><strong>Documentation »</strong></a>
    <br />
  </p>
</div>



<!-- TABLE OF CONTENTS 
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>
-->


<!-- ABOUT THE PROJECT -->
## About

About `wakeflow`.

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Installation

`wakeflow` may be most easily installed from the Python Package Index ([PyPI](https://example.com)), or can also be installed from the [GitHub repository](https://github.com/TomHilder/wakeflow) if you wish make contributions. Dependencies for `wakeflow` consist mostly of standard python libraries. We recommend using a package manager such as [Anaconda](https://www.anaconda.com/products/distribution) to make your life easier, but this is not required.



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



### PyPI (pip)

The easiest way to install `wakeflow` is via [PyPI](https://example.com), using `pip`:
```sh
pip install wakeflow
```
that's it!



### From source (GitHub)

If you want to contribute to, or modify `wakeflow`, you should install it from the [GitHub repository](https://github.com/TomHilder/wakeflow). Simply fork the repo using the button in the top right, and then clone it:
```sh
git clone https://github.com/<replace-by-your-username>/wakeflow.git
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

Use this space to show useful examples of how a project can be used. Additional screenshots, code examples and demos work well in this space. You may also link to more resources.

_For more examples, please refer to the [Documentation](https://example.com)_

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- CONTRIBUTING -->
## Contributing

Contributions to `wakeflow` are welcome. If you would like to implement a new feature, please:

1. Install using the above installatio from source instructions
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
[pypi-shield]: https://img.shields.io/badge/PyPI-0.0.1-brightgreen?style=for-the-badge&logo=appveyor
[pypi-url]: https://example.com
[python-shield]: https://img.shields.io/badge/Python-3.8%2B-orange?style=for-the-badge&logo=appveyor
[python-url]: https://www.python.org
[license-shield]: https://img.shields.io/badge/License-MIT-blue?style=for-the-badge&logo=appveyor
[license-url]: https://opensource.org/licenses/MIT
[JOSS-shield]: https://img.shields.io/badge/JOSS-coming%20soon-green?style=for-the-badge&logo=appveyor
[JOSS-url]: https://example.com

<!--
[product-screenshot]: images/screenshot.png
[Next.js]: https://img.shields.io/badge/next.js-000000?style=for-the-badge&logo=nextdotjs&logoColor=white
[Next-url]: https://nextjs.org/
[React.js]: https://img.shields.io/badge/React-20232A?style=for-the-badge&logo=react&logoColor=61DAFB
[React-url]: https://reactjs.org/
[Vue.js]: https://img.shields.io/badge/Vue.js-35495E?style=for-the-badge&logo=vuedotjs&logoColor=4FC08D
[Vue-url]: https://vuejs.org/
[Angular.io]: https://img.shields.io/badge/Angular-DD0031?style=for-the-badge&logo=angular&logoColor=white
[Angular-url]: https://angular.io/
[Svelte.dev]: https://img.shields.io/badge/Svelte-4A4A55?style=for-the-badge&logo=svelte&logoColor=FF3E00
[Svelte-url]: https://svelte.dev/
[Laravel.com]: https://img.shields.io/badge/Laravel-FF2D20?style=for-the-badge&logo=laravel&logoColor=white
[Laravel-url]: https://laravel.com
[Bootstrap.com]: https://img.shields.io/badge/Bootstrap-563D7C?style=for-the-badge&logo=bootstrap&logoColor=white
[Bootstrap-url]: https://getbootstrap.com
[JQuery.com]: https://img.shields.io/badge/jQuery-0769AD?style=for-the-badge&logo=jquery&logoColor=white
[JQuery-url]: https://jquery.com 
-->