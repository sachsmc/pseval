## Test environments
* local OS X install, R 3.2.3
* ubuntu 12.04 (on travis-ci), R 3.2.3
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There were 2 NOTEs:

* checking package dependencies ... NOTE
Package suggested but not available for checking: 'printr'

> printr is loaded if available in the vignette using "require". It makes table and help file output look pretty in the vignette, but it is optional. It is available from http://yihui.name/xran, which is listed in the Additional_repositories field of DESCRIPTION. 

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Michael C Sachs <sachsmc@gmail.com>’
Components with restrictions and base license permitting such:
  MIT + file LICENSE
File 'LICENSE':
  YEAR: 2015
  COPYRIGHT HOLDER: Michael C Sachs



## Downstream dependencies
There are no downstream dependencies.