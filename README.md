---
output:
  pdf_document: default
  html_document: default
  word_document: default
---
# YPmodelPhreg

This package fits a semi-parametric model to survival data with covariates. It provides covariate-adjusted comparison of two groups, referred to as the treatment group and the control group, respectively, of right censored data, where the binary group variable has separate short-term and long-term effects on the hazard function while effects of covariates, such as age and blood pressure, are proportional on the hazard. The model was studied in Yang and Prentice (2015), and it extends the two sample version of the short-term and long-term hazard ratio model proposed in Yang and Prentice (2005). The model extends the usual Cox proportional hazards model to allow more flexible hazard ratio patterns.

Let $T$ be the time-to-event outcome of interest, $Z$ the indicator that the observation is in the treatment group, and $W$ the covariate vector of dimension $p$. Let $\boldsymbol{U} = (\boldsymbol{Z},\boldsymbol{W}^\prime)^\prime$, where `$^\prime$' indicates transpose of a vector. The model assumes that the hazard function of $T$ given $\boldsymbol{U}$ is
\begin{align}
\lambda_{\boldsymbol{U}}(t) = \lambda_0(t)\frac{\exp(\boldsymbol{\alpha}^\prime \boldsymbol{W})}{\exp(-\beta_1 Z)S_0(t) + \exp(-\beta_2 Z)[1 - S_0(t)]},
\end{align}
where $\beta_1$ and $\beta_2$ are scalar parameters, $\boldsymbol{\alpha}$ is a $p$-dimensional parameter, and $\lambda_0(t)$ is an unknown baseline hazard function with corresponding survival function $S_0(t) =\exp\left[-\int_0^t\lambda_0(x)dx\right]$.

Under this model, the hazard ratio of the treatment group $(Z = 1)$ over the control group $(Z = 0)$ at covariate level $\boldsymbol{W} = \boldsymbol{w}$ is
\begin{align}
\frac{\lambda_{Z = 1, \boldsymbol{W} = \boldsymbol{w}}(t)}{\lambda_{Z = 0, \boldsymbol{W} = \boldsymbol{w}}(t)} = \frac{1}{\exp(-\beta_1 Z)S_0(t) + \exp(-\beta_2 Z)[1 - S_0(t)]}.
\end{align}

Similarly to the Cox proportional hazards model, the hazard ratio is free of the covariate $\boldsymbol{W}$. However, here the hazard ratio may be time-dependent. When $\beta_1 = \beta_2$, this model reduces to the Cox proportional hazards model. In general, the hazard ratio converges to the values $\exp(\beta_1)$ and $\exp(\beta_2)$, respectively, as $t$ tends to zero or the upper end of the support set of $S_0(t)$. Thus $\exp(\beta_1)$ and $\exp(\beta_2)$ can be interpreted as short-term and long-term hazard ratios, respectively. In addition to the proportional hazards pattern, other patterns of the hazard ratio can also be realized, including no initial effect, diminishing effect, crossing hazard or survival functions. Note that the model implies that the hazard ratio is monotone.

Let $h(t)$ be the hazard ratio function in (2). For a $\tau > 0$, the average hazard ratio over $[0,\tau]$ is
\begin{align}
AHR(\tau) = \int_0^\tau h(t)dt,
\end{align}
where the value of $\tau$ is controlled by the `tau` argument. The default value of $\tau$ is $0.9$ $\times$ the maximum of all observations, censored or not. However, the user can input a different value that is less than or equal to the maximum of all observations.

For a given data set, the package provides

* plots of estimated $h(t)$ and point-wise and simultaneous confidence bands of $h(t)$,
* point estimates and confidence intervals for the model parameters, $(\beta_1, \beta_2, \boldsymbol{\alpha}^\prime)^\prime$, and
* point estimate and confidence interval of $AHR(\tau)$. 

The default value of the confidence level is $95\%$, which can be controlled by the `alpha` argument (`alpha = 0.05` by default). The range $[L, U]$ of $t$ over which the simultaneous confidence band of h(t) is constructed can be chosen by the user. The default range $\left[L, U\right]$ is the time interval over which certain weight function is not too extreme (see Yang and Prentice, 2015 for details). When the user inputs a range $[L', U']$, then the package produces simultaneous confidence bands over the intersection of the default and user input intervals. The time points at which hazard ratios will be estimated along with the confidence intervals can be supplied by users via the ``time.hr`` argument. These time points must fall within the minimum and maximum values of uncensored observations.



## Installation
``` r
install.packages("YPmodelPhreg")
```
## Implementation
The main function in the package is the `ypreg` function. The `data` argument of the function must be a numeric matrix and consist of the variables for **time**, **event**, **group**, and **covariates** in this order:

* **time**: time until event or censoring
* **event**: censoring status (1 = event, 0 = censored)
* **group**: group indicator (1 = treatment, 0 = control)
* **covariates**: a set of numeric vectors of covariates

The additional arguments to `ypreg` are

* ``alpha``: a numeric value for the significance level for point-wise and simultaneous confidence bands (the default is $0.05$),
* ``time.hr``: a numeric vector of times at which hazard ratios will be estimated,
* ``L``: a numeric value for the lower bound of the range of $t$ over which the simultaneous confidence band of $h(t)$ is constructed,
* ``U``: a numeric value for the upper bound of the range of $t$ over which the simultaneous confidence band of $h(t)$ is constructed,
* ``repnum``: a numeric value for the number of replications for the re-sampling method, and
* ``tau``: a numeric value for the upper bound of $[0,\tau]$, over which the average hazard ratio will be obtained. 


## Example
The data set `colonexample` in the package is used as an example.

```r
library(YPmodelPhreg)
data(colonexample)
head(colonexample)

# fitting the model
res <- ypreg(colonexample, time.hr = c(1, 7))
res

# plotting the estimated hazard ratios with the confidence intervals
plot(res)
```

## Reference
Yang, S., & Prentice, R. (2005). Semiparametric analysis of short-term and long-term hazard ratios with two-sample survival data. *Biometrika*, 92(1), 1-17.

Yang, S., & Prentice, R. L. (2015). Assessing potentially time‐dependent treatment effect from clinical trials and observational studies for survival data, with applications to the Women's Health Initiative combined hormone therapy trial. *Statistics in medicine*, 34(11), 1801-1817.
