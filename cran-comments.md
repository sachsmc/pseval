This is a update that contains some important bugfixes and enhancements suggested by a review of a vignette submission to the R Journal.

## Test environments
* local Windows install, R 3.3.1
* ubuntu 12.04 (on travis-ci), R 3.3.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs. There was 1 WARNING, on R-devel:

* checking files in 'vignettes' ... WARNING
Files in the 'vignettes' directory newer than all files in 'inst/doc':
  'introduction.Rmd', 'psreferences.bib'
  
I don't understand this warning, as the files in inst/doc were created after the files in vignettes. I also see this warning on the CRAN R-devel checks but can't explain it or fix it. Sorry!

There was 1 NOTE:

Package suggested but not available for checking: 'printr'

printr is suggested but not required. It enhances the vignette appearence but the vignette runs without it. 

## Downstream dependencies
There are no downstream dependencies.