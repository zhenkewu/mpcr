How to install?
--------------
```r
install.packages(devtools)
devtools::install_github("zhenkewu/mpcr")
```
Why should someone use `mpcr`?
------------------------------

How does it compare to other existing solutions?
------------------------------------------------
- First R implementation of intervention effect estimation under matched-pair cluster randomized design, using baseline covariates for effect estimation consistency and efficiency.

- Reference is [here](http://onlinelibrary.wiley.com/doi/10.1111/biom.12214/full).

What are the main functions?
----------------------------
- `mpcr` function

Instructions:
-------------

-Prepare your dataset

You need to specify treatment arm, cluster number, primary outcome, covariates.
Note that suppose there are N pairs, then the first pair of clusters are numbered
as 1 and N+1; second pair as 2, N+2, etc.


Note
-----
current version does not have Bayes estimates, please specify BAYES_STATUS as FALSE. The Bayesian estimation code is ready. I will update it soon.
