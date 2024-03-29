
<!-- README.md is generated from README.Rmd. Please edit that file -->

# statstring

<!-- badges: start -->
<!-- badges: end -->

The goal of statstring is to facilitate formatting outputs of
statistical tests when they’re pasted into R markdown files in APA
style. Currently, statstring can format results of ANOVA outputs
generated by `stats::aov()`, `summary(stats::aov)`,
`rstatix::get_anova_table()`, and `apaTables::apa.aov.table()`. To
format the outputs for R markdown, simply pass one of the previously
listed ANVOA objects to the function `format_anova_string()`.

The package can also handle outputs from an independent samples
*t*-test, generated by `stats::t.test()`. The formatted output includes
the mean difference and a 95% confidence interval (CI) on the mean
difference.

You can also use the package to extract and format *R*<sup>2</sup> for a
model produced by `stats::lm()` or `summary(stats::lm())` using the
function `format_r2()`.

Finally, you can use this package to format point estimates and 95% CI
around the point estimate using `format_confint()`.

## Installation

You can install the released version of statstring from Github with:
`devtools::install_github("silverer/statstring")`

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(statstring)
data("warpbreaks")
#Results of a t-test
t.res = stats::t.test(breaks~wool, data = warpbreaks, var.equal=T)
t.res
#> 
#>  Two Sample t-test
#> 
#> data:  breaks by wool
#> t = 1.6335, df = 52, p-value = 0.1084
#> alternative hypothesis: true difference in means between group A and group B is not equal to 0
#> 95 percent confidence interval:
#>  -1.319679 12.875235
#> sample estimates:
#> mean in group A mean in group B 
#>        31.03704        25.25926
t.str = format_tstat_apa(t.res)
#Results of ANOVA
aov.res = stats::aov(breaks~wool*tension, data = warpbreaks)
summary(aov.res)
#>              Df Sum Sq Mean Sq F value   Pr(>F)    
#> wool          1    451   450.7   3.765 0.058213 .  
#> tension       2   2034  1017.1   8.498 0.000693 ***
#> wool:tension  2   1003   501.4   4.189 0.021044 *  
#> Residuals    48   5745   119.7                     
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#Getting all F-statistics
all.aov.results = format_anova_string(summary(aov.res))

#Results of linear model
lm.res = stats::lm(breaks~ wool + tension, data = warpbreaks)
summary(lm.res)
#> 
#> Call:
#> stats::lm(formula = breaks ~ wool + tension, data = warpbreaks)
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -19.500  -8.083  -2.139   6.472  30.722 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)   39.278      3.162  12.423  < 2e-16 ***
#> woolB         -5.778      3.162  -1.827 0.073614 .  
#> tensionM     -10.000      3.872  -2.582 0.012787 *  
#> tensionH     -14.722      3.872  -3.802 0.000391 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 11.62 on 50 degrees of freedom
#> Multiple R-squared:  0.2691, Adjusted R-squared:  0.2253 
#> F-statistic: 6.138 on 3 and 50 DF,  p-value: 0.00123
r2 = format_r2(lm.res)
```

For the t-test, the output looks like *t*(52) = 1.63, *p* = .11,
M<sub>diff</sub> (95% CI) = 5.78 (-1.32, 12.88).

For ANOVAs, the output corresponds to the order in which terms were
entered into the model. For example, the term for wool (the first factor
in the model) is: *F*(1, 48) = 3.77, *p* = .06, for tension is *F*(2,
48) = 8.50, *p* \< .001, and for the interaction between tension and
wool: *F*(2, 48) = 4.19, *p* = .02

For a linear model using breaks as the criterion and wool and tension as
predictors, *R*<sup>2</sup> = 0.27

You can also get formatted text outputs that do not include the special
characters used for formatting things in markdown, such as the
underscore for italics. To turn off markdown formatting, pass the
argument `as.markdown=F` to your function.
