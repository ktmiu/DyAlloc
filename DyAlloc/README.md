
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DyAlloc

<!-- badges: start -->
<!-- badges: end -->

The goal of DyAlloc is to provide clinical trial group allocation for a
new subject entering the study. Imbalance scores for each possible
treatment group for the new subject are calculated given all the
profiles from the previous subjects, along with the treatment allocation
probability for each group.

## Installation

You can install the released version of DyAlloc from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("DyAlloc")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ktmiu/DyAlloc")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(DyAlloc)
## basic example code
#generate a data frame of 50 subjects in a 4 treatment group study
df1 <- data.frame("gender" = sample(c("female", "male"),50, TRUE),
                 "age" = sample(c("0-35", ">35"),50, TRUE),
                 "race" = sample(c("black", "white"),50,TRUE),
                 "trauma severity"= sample(c("good", "bad"),50,TRUE),
                 "treatment"=sample(c("1","2","3","4"),50, TRUE,c(1,1,5,10)),
                  stringsAsFactors = TRUE)
#covariates of the study
covars <-c("gender","age","race","trauma severity")
#new subject's profile
newsub <-c("male","0-35","black","bad")
#weight of the overall study, the stratum, and the 4 factors 
weight <-c(1,1,rep(1,4))
#print out the final weighted imbalance scores
imb <-WeiImb(df1,4,covars,weight,newsub,FALSE)
imb
#> [1] 6.954902 6.454902 6.464706 7.127451
#probabilities for the new subject to be assigned to the each of the treatment group 
Trt_Prob(imb,0.2)
#> [1] 0.0 0.8 0.2 0.0

#Using CAR function to get the final treatment assignment
DyAlloc(df1,4,covars,weight,newsub,FALSE)
#> [1] "2"
```
