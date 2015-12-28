[![Travis-CI Build Status](https://travis-ci.org/sachsmc/pseval.svg?branch=master)](https://travis-ci.org/sachsmc/pseval)

# pseval: Methods for Evaluting Principal Surrogates of Treatment Response

## Installation

`pseval` is an R package aimed at implementing existing methods for surrogate evaluation using a flexible and common interface. Development will take place  on [the Github page](https://github.com/sachsmc/pseval), and the current version of the package can be installed as shown below. First you must install the `devtools` package, if you haven't already `install.packages("devtools")`. 

```{r eval = FALSE}
devtools::install_github("sachsmc/pseval")
```

Check out the [vignette](https://sachsmc.github.io/pseval) for methodological details and information on how to use the package.

## References


- [Gabriel and Gilbert, 2014. _Evaluating principal surrogate endpoints with time-to-event data accounting for time-varying treatment efficacy_](http://biostatistics.oxfordjournals.org/content/15/2/251)
- [Huang and Gilbert, 2011. _Comparing Biomarkers as Principal Surrogate Endpoints_](http://onlinelibrary.wiley.com/doi/10.1111/j.1541-0420.2011.01603.x/full)
- [Gilbert and Hudgens, 2008. _Evaluating Candidate Principal Surrogate Endpoints_](http://onlinelibrary.wiley.com/doi/10.1111/j.1541-0420.2008.01014.x/full)
- [Huang, Gilbert, and Wolfson, 2013. _Design and Estimation for Evaluating Principal Surrogate Markers in Vaccine Trials_](http://onlinelibrary.wiley.com/doi/10.1111/biom.12014/full)