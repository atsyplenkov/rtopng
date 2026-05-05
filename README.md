# rtopng

<!-- badges: start -->
[![R-CMD-check](https://github.com/atsyplenkov/rtopng/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/atsyplenkov/rtopng/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`rtopng` is a fork and continuation of the original [`rtop` package](https://CRAN.R-project.org/package=rtop). It keeps the core top-kriging functionality for interpolation with variable spatial support while modernising the package internals and extending the feature set.

This fork is intended to be rewritten using more modern R tooling, including `mirai` for parallel processing, `S7` for a more explicit object system, and `stars` for spatiotemporal datasets in place of the deprecated `spacetime` package. It also aims to add broader support for universal kriging (hope so!).

## Installation

You can install the development version of `rtopng` from GitHub with:

``` r
# install.packages("pak")
pak::pak("atsyplenkov/rtopng")
```

## Roadmap

- [x] Migrate tests to testthat
- [ ] Replace spacetime with stars
- [ ] Add a proper OOP structure through S7
- [ ] Remove data.table and reshape2
- [ ] Add mirai for parallel processing
- [ ] Add universal kriging support

## Acknowledgements

`rtopng` builds on the original `rtop` package by Jon Olav Skøien and the method description published in Skøien et al. ([2014](doi:10.1016/j.cageo.2014.02.009)).

## License

`rtopng` is distributed under GPL-3, consistent with the `rtop`. For the upstream package lineage and original package reference, see the original [`rtop` CRAN page](https://CRAN.R-project.org/package=rtop).