<!-- badges: start -->
[![Travis build status](https://travis-ci.org/ntthung/ldsr.svg?branch=master)](https://travis-ci.org/ntthung/ldsr)
[![GPL license](https://img.shields.io/badge/License-GPL-blue.svg)](http://perso.crans.org/besson/LICENSE.html)
[![GitHub version](https://badge.fury.io/gh/ntthung%2Fldsr.svg)](https://badge.fury.io/gh/ntthung%2Fldsr)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/ldsr)](https://cran.r-project.org/package=ldsr)
<!-- badges: end -->


## Streamflow Reconstruction Using Linear Dynamical System

A typical streamflow reconstruction model is linear: the relationship between streamflow *y* and the paleoclimate proxies *u* is modelled as

<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=y_t = \alpha %2B \beta u_t %2B \varepsilon_t">
</p>

Linear models do not account for the catchment state and the catchment memory effect. To model these, we experimented with using the Linear Dynamical System (LDS) model. We model the relationship between streamflow and paleoclimate proxies as follows:

<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign*%7D%20x_%7Bt%2B1%7D%20%3D%20Ax_t%20%2B%20Bu_t%20%2B%20w_t%5C%5C%20y_t%20%3D%20Cx_t%20%2B%20Du_t%20%2B%20v_t%20%5Cend%7Balign*%7D%20">
</p>

Observe that linear regression is a special case of LDS. The constant term ![\alpha](https://render.githubusercontent.com/render/math?math=%5Calpha) of linear regression is replaced by a state-dependent term ![Cx_t](https://render.githubusercontent.com/render/math?math=Cx_t), and the system state follows a state transition equation.

We described the method in full detail, together with a case study for the Ping River (Thailand), in Nguyen and Galelli (2018)---please cite this paper when you use the package. We also used the method to reconstruct streamflow for 48 stations in 16 countries in Asia (Nguyen et al, 2019). 

For a tutorial on how to use the package, please see the [package's vignette](https://cran.r-project.org/web/packages/ldsr/vignettes/ldsr.html)

## Installation

LDS is somewhat computationally heavy. To speed it up, I wrote the core routines in C++, which needs compilation. On Windows, you need to first install Rtools. [For R 4.0.0](https://cran.r-project.org/bin/windows/Rtools/). [For older R](https://cran.r-project.org/bin/windows/Rtools/history.html).

`ldsr` is now available on CRAN, so you can install with

```
install.packages('ldsr')
```

To install the development version

```
install.packages('remotes')
remotes::install_github('ntthung/ldsr', build_vignettes = TRUE)
```

## References

Nguyen, H. T. T., & Galelli, S. (2018). A linear dynamical systems approach to streamflow reconstruction reveals history of regime shifts in northern Thailand. Water Resources Research, 54, 2057– 2077. https://doi.org/10.1002/2017WR022114 

Nguyen, H. T. T., Turner, S. W., Buckley, B. M., & Galelli, S. (2019). Coherent streamflow variability in Monsoon Asia over the past eight centuries---links to oceanic drivers. https://doi.org/10.31223/osf.io/5tg68
