# poisson.mash.alpha

An R package for Poisson multivariate adaptive shrinkage ("Poisson
mash"). The Poisson mash model and methods are described in
[Overleaf][overleaf].

## Quick Start

Install and load the package using [remotes][remotes]:

```R
install.packages("remotes")
remotes::install_github("stephenslab/poisson.mash.alpha")
library(poisson.mash.alpha)
```

Note that installing the package will require a C++ compiler setup
that is appropriate for the version of R installed on your
computer. For details, refer to the documentation on the
[CRAN website][cran].

Note this should automatically install the ashr package; if it does
not, please go [here][ashr].

## License

Copyright (c) 2021-2022, Yusha Liu, Peter Carbonetto and Matthew
Stephens.

All source code and software in this repository are made available
under the terms of the [MIT license][mit-license].

[mit-license]: https://opensource.org/licenses/mit-license.html
[remotes]: https://github.com/r-lib/remotes
[overleaf]: https://www.overleaf.com/read/pqzccwzpqpvs
[ashr]: https://github.com/stephens999/ashr
[cran]: https://cran.r-project.org
