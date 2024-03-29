# GRPtests: Goodness-of-fit testing for high-dimensional generalized linear models

## Contents of this repository

This repository contains:

1. R package GRPtests in the [GRPtests_Rpackage_update_Jan_2021](https://github.com/jankova/GRPtests/tree/master/GRPtests_Rpackage_update_Jan_2021) folder. 
The package can also be installed from CRAN (see the Installation section below).

2. Code reproducing empirical results from Section 4 of [[1]](https://rss.onlinelibrary.wiley.com/doi/10.1111/rssb.12371).

## Installation of R package GRPtests

Installation in R from CRAN repository:

```
install.packages("GRPtests")

library("GRPtests")
```

Installation in R from this github repository:

```
install.packages(devtools)

library(devtools)

install_github("jankova/GRPtests/GRPtests_Rpackage_update_Jan_2021/GRPtests")
```


## Examples from empirical section of [1]

The code reproducing empirical results from Section 4 of [1] is available in <br/><br/>
[Example_Section_4-1](https://github.com/jankova/GRPtests/blob/master/Example_Section_4-1.R),<br/>
[Example_Section_4-2](https://github.com/jankova/GRPtests/blob/master/Example_Section_4-1.R),<br/>
[Example_Section_4-3](https://github.com/jankova/GRPtests/blob/master/Example_Section_4-1.R),<br/>
[Example_Section_4-4](https://github.com/jankova/GRPtests/blob/master/Example_Section_4-1.R).<br/>

## Method heuristics
We briefly sketch the method's heuristics. 

Let Y be the target vector and let X be the matrix with features as columns. 

Split the observations randomly into two parts (X_A, y_A) and (X_B, y_B). 

Fit a GLM regression of y_A on X_A and y_B on X_B and for both fitted regressions, compute the Pearson residuals
R_A and R_B.

The main idea of the method is then as follows: if the logistic regression model was not a good fit, we would expect that some nonlinear signal was left in
the residuals.
Therefore we use a ML method (by default the random forest, but any method may be used) to predict the leftover signal from residuals R_A.

Using the random forest to fit R_A on X_A, we obtain a prediction function f_A(). 
If there was indeed nonlinear signal left in the residuals and the random forest picked it up,
then this signal should also be present in the residuals from part B. 
Thus if we predict the random forest f_A() on X_B and compute the scalar product of f_A(X_B) and R_B, this would be large and we would reject the null hypothesis.

A schematic illustration of the procedure is pictured below.

<img src="grpimage.jpg" alt="methodology_diagram" width="300"/>

## References

[1] [Janková, J., Shah, R. D., Bühlmann, P. and Samworth, R. J., Goodness-of-fit testing in high-dimensional generalized linear models (2020), Journal of the Royal Statistical Society 82, Part 3, pp. 773–795](https://rss.onlinelibrary.wiley.com/doi/10.1111/rssb.12371)
