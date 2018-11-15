# ldsr
Streamflow reconstruction using linear dynamical system

To install this package, type the following commands into your RStudio console:

```
install.packages('devtools')
devtools::install_github('ntthung/ldsr', build_opts = c("--no-resave-data", '--no-manual'))
```

Note that if you don't already have Rtools (a set of tools to build packages from sources), RStudio will prompt you to install it. Just follow the instructions. `ldsr` will be installed automatically after Rtools.

To view the package's vignette, type

`browseVignettes('ldsr')`
