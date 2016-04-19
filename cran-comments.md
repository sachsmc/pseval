This is a resubmission. The following change was made: 

Bug fix to deal with error in vignette. 

## Test environments
* local OS X install, R 3.2.4
* ubuntu 12.04 (on travis-ci), R 3.2.4
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There were 2 NOTEs:

* checking package dependencies ... NOTE
Package suggested but not available for checking: 'printr'

> printr is loaded conditionally in the vignette using "require". It makes some of the output look cleaner in the vignette, but it is optional and will build without printr being available. 

Possibly mis-spelled words in DESCRIPTION:
  counterfactual (14:38)

> This is not a mis-spelling. It is a commonly used term in the causal inference field. 

## Downstream dependencies
There are no downstream dependencies.