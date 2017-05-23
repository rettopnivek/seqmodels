# seqmodels

A R package with density, distribution, quantile, and random number generation functions for a collection of sequential sampling and response time models.

## Getting started

### Prerequisites

The program R ( version >= 3.0 )
The R package ['Rcpp'](https://cran.r-project.org/web/packages/Rcpp/index.html)
The R package ['RcppParallel'](https://rcppcore.github.io/RcppParallel/)

### Installation

To easily install the 'seqmodels' package, you'll need the 'devtools' package:  
```
install.packages("devtools")
library(devtools)
```

The 'seqmodels' package can then be installed via the following command:  
```
install_github("rettopnivek/seqmodels")
```

## Using the package

To load the package:
```
library(seqmodels)
```

To list the available functions:
```
ls(pos = "package:seqmodels")
```

Details and examples for a specific function can be obtained via:
```
help( "function" )
```
subsituting the appropriate function name for "function".

## Authors

Kevin Potter

## License

MIT
