This is a resubmission. The following change was made: 

https://cran.r-project.org/web/packages/survey was changed to https://cran.r-project.org/package=survey in inst/doc/introduction.html.

## Test environments
* local OS X install, R 3.2.3
* ubuntu 12.04 (on travis-ci), R 3.2.3
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking package dependencies ... NOTE
Package suggested but not available for checking: 'printr'

> printr is loaded conditionally in the vignette using "require". It makes some of the output look cleaner in the vignette, but it is optional and will build without printr being available. 

## Downstream dependencies
There are no downstream dependencies.