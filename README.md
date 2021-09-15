# poisson.mash.alpha

An R package for Poisson multivariate adaptive shrinkage ("Poisson
mash"). The Poisson mash model and methods are described in
[Overleaf][overleaf].

## Quick Start

Clone or download the git repository from GitHub, then install the
latest version of the package using [remotes][remotes]:

```R
install.packages("remotes")
remotes::install_local("poisson.mash.alpha")
```

This assumes that your R working directory contains a local copy of
the git repository. This should automatically install the ashr
package; if it does not, please go [here][ashr].

## License

Copyright (c) 2021, Yusha Liu, Peter Carbonetto and Matthew
Stephens.

All source code and software in this repository are made available
under the terms of the [MIT license][mit-license].

[mit-license]: https://opensource.org/licenses/mit-license.html
[remotes]: https://github.com/r-lib/remotes
[overleaf]: https://www.overleaf.com/read/pqzccwzpqpvs
[ashr]: https://github.com/stephens999/ashr
