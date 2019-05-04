
# Bayesian Statistics Final Project Report

## Efe Atikkan

The aim of this project is to examine the effects of seed type and root
extract on germination. Bayesian analysis will be done for the given
data and a suitable model will be fit. Diagnostics about markov Chain
Monte Carlo and comparison with some other possible models will be
presented.

# The Data

The data used in this project consist of porportion of seeds that are
germinated in 21 plates. The data set is arranged in a 2 x 2 factor
layout and these factors are: - Type of the seed : O. aegyptiaca 75 and
O. aegyptiaca 73 - Type of the root extract: Bean or Cucumber  
Total number of seeds (n), number of germinated seeds (r) and their
ratio is given for each plate .

seed\_type=x1 seed\_extract=x2

``` r
"N" <- 21
"r" <-
  c(10, 23, 23, 26, 17, 5, 53, 55, 32, 46, 10, 8, 10, 8, 23, 0, 
    3, 22, 15, 32, 3)
"n" <-
  c(39, 62, 81, 51, 39, 6, 74, 72, 51, 79, 13, 16, 30, 28, 45, 
    4, 12, 41, 30, 51, 7)
"x1" <-
  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1)
"x2" <-
  c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 
    1)
```

# The model

A random effect logistic model is considered for this data.

\[r_i \sim Binomial (p_i, n_i)\]

where $ r\_i $ is the number of germinated seeds in plate i, $ n\_i $ is
the total number of seeds in plate i and p\_i is the probability of
germination on plate i. and $ p\_i $ is defined
as:

\[logit(p_i) = \alpha_0 + \alpha_1*x_{1i} + \alpha_2 * x_{2i} + \alpha_{12} * x_{1i}*x_{2i} + b_i \]

Here b\_i is the random effect distributed as:
\[b_i \sim normal(0, \tau)\]

PRIORS

there are 4 variables that we need to give prior distributions:
\[\alpha_0, \alpha_1, \alpha_2, \alpha_{12} , \tau\]

Since he don’t any more information about this dataset and this topic,
it is better to give non-informative priors. Normal distribution with
mean 0 and very high standard deviation would be suitable for alpha
variables.

\[\alpha_0 \sim dnorm(0, 100^2)\] \(\alpha_1 \sim dnorm(0, 100^2)\)
\(\alpha_2 \sim dnorm(0, 100^2)\) \(\alpha_12 \sim dnorm(0, 100^2)\)
\(\tau \sim dgamma(0.1, 0.1)\)

plotting of priors:

``` r
curve(dnorm(x,0,1000),from=-5000,5000,main='alpha prior dist.')
```

![](project_markdown_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
curve(dgamma(x,0.01,0.01))
```

![](project_markdown_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

Initialization of parameters:

``` r
ino=list(list("tau" = 1, "alpha0" =
  0,
"alpha1" =
  0,
"alpha2" =
  0,
"alpha12" =
  0))

ino2=list(list("tau" = 0.5, "alpha0" =
  1,
"alpha1" =
  1,
"alpha2" =
  1,
"alpha12" =
  1))

ino3=list(list("tau" = 0.1, "alpha0" =
  -1,
"alpha1" =
  -1,
"alpha2" =
  -1,
"alpha12" =
  -1))
```

# Markov Chain Monte Carlo Using BUGS

The model explained in previous parts will be implemented using
R2WinBUGS package. The models are saved to a txt file for further usage.

Note: BUGS uses tau (precision ) parameter instead of variancec
(sigma^2) which defined as tau=1/(sigma^2).

The parameters of the BUGS model are defined as following:

``` r
my.data = list(N = N, n = n, r = r, x1= x1 , x2= x2)
params = c("alpha0", "alpha1", "alpha2","alpha12", "sigma")
```

BUGS Model is
runned:

``` r
bugsModel1 = bugs( model.file="m1.bugs" , data=my.data , inits=c(ino,ino2,ino3) , 
                        n.chains=3,parameters.to.save = params ,
                        n.iter =20000,n.thin=1,n.burnin = 1000,bugs.directory = "C://Users//atief//Desktop//efe//master//derslerr//2nd semester//bayesian//winbugs//WinBUGS14")
```

Results Summary:

``` r
print(bugsModel1)
```

    ## Inference for Bugs model at "m1.bugs", fit using WinBUGS,
    ##  3 chains, each with 20000 iterations (first 1000 discarded)
    ##  n.sims = 57000 iterations saved
    ##           mean  sd 2.5%  25%   50%   75% 97.5% Rhat n.eff
    ## alpha0    -0.6 0.2 -1.0 -0.7  -0.6  -0.4  -0.1    1  2800
    ## alpha1     0.1 0.3 -0.6 -0.1   0.1   0.3   0.7    1  8200
    ## alpha2     1.4 0.3  0.8  1.2   1.4   1.5   1.9    1  2200
    ## alpha12   -0.8 0.5 -1.8 -1.1  -0.8  -0.5   0.1    1  5000
    ## sigma      0.3 0.1  0.1  0.2   0.3   0.4   0.6    1  4000
    ## deviance 100.5 6.3 89.6 96.0 100.1 104.6 113.7    1 15000
    ## 
    ## For each parameter, n.eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
    ## 
    ## DIC info (using the rule, pD = Dbar-Dhat)
    ## pD = 11.7 and DIC = 112.2
    ## DIC is an estimate of expected predictive error (lower deviance is better).

``` r
bugsModel1$mean
```

    ## $alpha0
    ## [1] -0.5512509
    ## 
    ## $alpha1
    ## [1] 0.07645881
    ## 
    ## $alpha2
    ## [1] 1.360139
    ## 
    ## $alpha12
    ## [1] -0.8402492
    ## 
    ## $sigma
    ## [1] 0.3229239
    ## 
    ## $deviance
    ## [1] 100.5292

``` r
for (i in 1:5){
plot(bugsModel1$sims.array[,1,i], type = "l", main=names(bugsModel1$sims.list)[i], col='red')
  
}
```

![](project_markdown_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-10-3.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-10-4.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-10-5.png)<!-- -->

There shouldn’t be correlation between the samples of parameter
distributions. Autocorrelation graphs help us to check this situation
and gives us the visual inspection of correlation versus lagged sample
values.

Lets look at the autocorrelation graphs

``` r
par(mfrow=c(3,2))
for (i in 1:5){
  acf(bugsModel1$sims.array[,1,i],main=names(bugsModel1$sims.list)[i])
}
par(mfrow=c(1,1))
```

![](project_markdown_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

We can see that the correlation between samples decreases in a pretty
fast manner. It is another check to make sure our MCMC model converged
well and no sign of problematic issue which may effect the inference.

# Density Plots

``` r
for (i in 1:5 ){
plot(density(bugsModel1$sims.array[,1,i]),main =paste("Posterior Density Plot of ", names(bugsModel1$sims.list)[i]), lwd=2)
}
```

![](project_markdown_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-12-3.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-12-4.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-12-5.png)<!-- -->

# Frequentist logistic reg.

We can estimate the parameters of the model using frequentist approach.

    ## 
    ## Call:  glm(formula = vv ~ x1 + x1 * x2 + x2, family = binomial(link = logit))
    ## 
    ## Coefficients:
    ## (Intercept)           x1           x2        x1:x2  
    ##     -0.5262      -0.2000       1.4479      -0.8478  
    ## 
    ## Degrees of Freedom: 20 Total (i.e. Null);  17 Residual
    ## Null Deviance:       3.911 
    ## Residual Deviance: 1.815     AIC: 29.94

``` r
bugsModel1$mean
```

    ## $alpha0
    ## [1] -0.5512509
    ## 
    ## $alpha1
    ## [1] 0.07645881
    ## 
    ## $alpha2
    ## [1] 1.360139
    ## 
    ## $alpha12
    ## [1] -0.8402492
    ## 
    ## $sigma
    ## [1] 0.3229239
    ## 
    ## $deviance
    ## [1] 100.5292

As we can see, parameter estimates inferred by frequentist logistic
regression are very close to the means of the distributions estimated by
bayesian model.

High Density Intervals for posterior distributions are calculated.

Given a credibility mass, HDI values represents the interval which
highest density is located in the given distribution.

``` r
library(HDInterval)
hdi_1=hdi(as.mcmc(bugsModel1$sims.array[,1,1]))
hdi_2=hdi(as.mcmc(bugsModel1$sims.array[,1,2]))
hdi_3=hdi(as.mcmc(bugsModel1$sims.array[,1,3]))
hdi_4=hdi(as.mcmc(bugsModel1$sims.array[,1,4]))
```

In the following graphs, posterior distributions of estimated parameters
are plotted. The blue lines represents the Highest Density Intervals and
the red line represents the frequentist point estimate of the
parameters.

``` r
for (i in 1:4){
plot(density(bugsModel1$sims.array[,1,i]),main =paste("Posterior Density Plot of ", names(bugsModel1$sims.list)[i]),lwd=2)
abline(v=f_log_reg$coefficients[i],col='red', lwd=2)
abline(v=eval(parse(text=sprintf("hdi_%d",i)))[1],col='blue', lwd=1.5)
abline(v=eval(parse(text=sprintf("hdi_%d",i)))[2],col='blue', lwd=1.5)
}
```

![](project_markdown_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-16-3.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-16-4.png)<!-- -->

On possible way of further developing our bayesian analysis is to
generate data from our base model and fit it to several other methods in
order to (hopefully) observe that it fits best to the model that
generated it. In this section we will go over these methodology step by
step.

# Data Generation from the base model

``` r
bugsModel1$mean
```

    ## $alpha0
    ## [1] -0.5512509
    ## 
    ## $alpha1
    ## [1] 0.07645881
    ## 
    ## $alpha2
    ## [1] 1.360139
    ## 
    ## $alpha12
    ## [1] -0.8402492
    ## 
    ## $sigma
    ## [1] 0.3229239
    ## 
    ## $deviance
    ## [1] 100.5292

``` r
alpha0_val=bugsModel1$mean$alpha0
alpha1_val=bugsModel1$mean$alpha1
alpha2_val=bugsModel1$mean$alpha2
alpha12_val=bugsModel1$mean$alpha12
sigma_val=bugsModel1$mean$sigma

x1_vals=c()
x2_vals=c()

n_vals= c()
k=c()
for (i in 1:200){
  x1_r=sample(c(0,1),1, replace = T)
  x2_r=sample(c(0,1),1, replace=T)
  n_val= sample(3:55, 1, replace = T)
  
  x1_vals= c(x1_vals,x1_r)
  x2_vals= c(x2_vals,x2_r)
  n_vals= c(n_vals, n_val)
  
  
  vv=alpha0_val + alpha1_val * x1_r + alpha2_val * x2_r + alpha12_val * x1_r * x2_r +rnorm(1, 0,sigma_val)
  #print(inv.logit(vv))
  k=c(k,rbinom(1, n_val, inv.logit(vv)))
    
}
```

``` r
generated_data_sample=cbind(as.data.frame(n_vals[1:10]),as.data.frame(k[1:10]),
                            as.data.frame(x1_vals[1:10]),as.data.frame(x2_vals[1:10])
                            )
names(generated_data_sample)=c('total no', 'germinated', 'x1', 'x2')
print(generated_data_sample)
```

    ##    total no germinated x1 x2
    ## 1        45         20  1  0
    ## 2        42         15  1  1
    ## 3        43         33  0  1
    ## 4        30         21  0  1
    ## 5        14          4  0  0
    ## 6        53         29  1  1
    ## 7        48         37  0  1
    ## 8        13         10  1  1
    ## 9        41         21  1  0
    ## 10       13          5  0  0

``` r
"N_F" <- 200
"r_F" <- k
"n_F" <- n_vals
"x1_F" <- x1_vals
"x2_F" <- x2_vals
```

# Creation and fitting of the first (base) model

``` r
ino=list(list("tau" = 1, "alpha0" =0,"alpha1" =0,"alpha2" =0,"alpha12" =0))
ino2=list(list("tau" = 0.5, "alpha0" =1,"alpha1" =1,"alpha2" =1,"alpha12" =1))
ino3=list(list("tau" = 0.1, "alpha0" =-1,"alpha1" =-1,"alpha2" =-1,"alpha12" =-1))


modelString = "
      model
      {
         for( i in 1 : N ) {
            r[i] ~ dbin(p[i],n[i])
            b[i] ~ dnorm(0.0,tau)
            logit(p[i]) <- alpha0 + alpha1 * x1[i] + alpha2 * x2[i]  + alpha12 * x1[i] * x2[i]  + b[i]
         }
         alpha0 ~ dnorm(0.0,1.0E-4)
         alpha1 ~ dnorm(0.0,1.0E-4)
         alpha2 ~ dnorm(0.0,1.0E-4)
         alpha12 ~ dnorm(0.0, 1.0E-4)
         tau ~ dgamma(0.01,0.01)
         sigma <- 1 / sqrt(tau)
      }"

writeLines( modelString , con="fake_m1.bugs" )
```

The parameters of the BUGS model are defined as following:

``` r
my.data = list(N = N_F, n = n_F, r = r_F, x1= x1_F , x2= x2_F)
params = c("alpha0", "alpha1", "alpha2", "alpha12" , "sigma")
```

Bugs Model is
runned:

``` r
fake_data_model1 = bugs( model.file="fake_m1.bugs" , data=my.data , inits=c(ino,ino2,ino3) , 
                        n.chains=3,parameters.to.save = params ,
                        n.iter =10000,n.thin=1,n.burnin = 2000,bugs.directory = "C://Users//atief//Desktop//efe//master//derslerr//2nd semester//bayesian//winbugs//WinBUGS14")
```

Model 1 summary

``` r
fake_data_model1
```

    ## Inference for Bugs model at "fake_m1.bugs", fit using WinBUGS,
    ##  3 chains, each with 10000 iterations (first 2000 discarded)
    ##  n.sims = 24000 iterations saved
    ##           mean   sd  2.5%   25%   50%   75% 97.5% Rhat n.eff
    ## alpha0    -0.6  0.1  -0.8  -0.7  -0.6  -0.6  -0.5    1  2300
    ## alpha1     0.1  0.1  -0.1   0.0   0.1   0.2   0.3    1  2600
    ## alpha2     1.4  0.1   1.2   1.3   1.4   1.5   1.6    1  1400
    ## alpha12   -0.7  0.2  -1.0  -0.8  -0.7  -0.6  -0.3    1  2800
    ## sigma      0.4  0.0   0.3   0.3   0.4   0.4   0.5    1  7500
    ## deviance 900.6 19.0 864.8 887.7 900.0 913.1 939.3    1  9000
    ## 
    ## For each parameter, n.eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
    ## 
    ## DIC info (using the rule, pD = Dbar-Dhat)
    ## pD = 92.6 and DIC = 993.2
    ## DIC is an estimate of expected predictive error (lower deviance is better).

# More basic model

In this model, a less complex relation is built. The cross-product term
(x1\*x2) alpha\_12 is not fitted.

``` r
ino=list(list("tau" = 1, "alpha0" =
  0,
"alpha1" =
  0,
"alpha2" =
  0))

ino2=list(list("tau" = 0.5, "alpha0" =1,"alpha1" =1,"alpha2" =1))
ino3=list(list("tau" = 0.1, "alpha0" =-1,"alpha1" =-1,"alpha2" =-1))


modelString = "
      model
      {
         for( i in 1 : N ) {
            r[i] ~ dbin(p[i],n[i])
            b[i] ~ dnorm(0.0,tau)
            logit(p[i]) <- alpha0 + alpha1 * x1[i] + alpha2 * x2[i] +b[i]
         }
         alpha0 ~ dnorm(0.0,1.0E-4)
         alpha1 ~ dnorm(0.0,1.0E-4)
         alpha2 ~ dnorm(0.0,1.0E-4)
         tau ~ dgamma(0.01,0.01)
         sigma <- 1 / sqrt(tau)
      }"

writeLines( modelString , con="fake_m2.bugs" )

my.data = list(N = N_F, n = n_F, r = r_F, x1= x1_F , x2= x2_F)
params = c("alpha0", "alpha1", "alpha2" , "sigma")

fake_data_model2 = bugs( model.file="fake_m2.bugs" , data=my.data , inits=c(ino,ino2,ino3) , 
                        n.chains=3,parameters.to.save = params ,
                        n.iter =10000,n.thin=1,n.burnin = 1000,bugs.directory = "C://Users//atief//Desktop//efe//master//derslerr//2nd semester//bayesian//winbugs//WinBUGS14")
```

``` r
fake_data_model2
```

    ## Inference for Bugs model at "fake_m2.bugs", fit using WinBUGS,
    ##  3 chains, each with 10000 iterations (first 1000 discarded)
    ##  n.sims = 27000 iterations saved
    ##           mean   sd  2.5%   25%   50%   75% 97.5% Rhat n.eff
    ## alpha0    -0.5  0.1  -0.6  -0.5  -0.5  -0.4  -0.3    1  1600
    ## alpha1    -0.2  0.1  -0.4  -0.3  -0.2  -0.2  -0.1    1  1100
    ## alpha2     1.0  0.1   0.9   1.0   1.0   1.1   1.2    1  2300
    ## sigma      0.4  0.0   0.3   0.4   0.4   0.4   0.5    1  7100
    ## deviance 904.5 19.4 868.2 891.0 903.9 917.2 944.1    1 23000
    ## 
    ## For each parameter, n.eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
    ## 
    ## DIC info (using the rule, pD = Dbar-Dhat)
    ## pD = 98.2 and DIC = 1002.7
    ## DIC is an estimate of expected predictive error (lower deviance is better).

As we can see from DIC part, the first base model (which the data is
generated from) has a lower DIC value, which means it’s better fit to
the data as expected.

# Least complex model

Another less complex model is fitted here; only x1 variable is used.

``` r
ino=list(list("tau" = 1, "alpha0" =0,"alpha1" =0))
ino2=list(list("tau" = 0.5, "alpha0" =1,"alpha1" =1))
ino3=list(list("tau" = 0.1, "alpha0" =-1,"alpha1" =-1))

modelString = "
      model
      {
         for( i in 1 : N ) {
            r[i] ~ dbin(p[i],n[i])
            b[i] ~ dnorm(0.0,tau)
            logit(p[i]) <- alpha0 + alpha1 * x1[i] +b[i]
         }
         alpha0 ~ dnorm(0.0,1.0E-4)
         alpha1 ~ dnorm(0.0,1.0E-4)
         tau ~ dgamma(0.01,0.01)
         sigma <- 1 / sqrt(tau)
      }"

writeLines( modelString , con="fake_m3.bugs" )

my.data = list(N = N_F, n = n_F, r = r_F, x1= x1_F)
params = c("alpha0", "alpha1" , "sigma")

fake_data_model3 = bugs( model.file="fake_m3.bugs" , data=my.data , inits=c(ino,ino2,ino3) , 
                        n.chains=3,parameters.to.save = params ,
                        n.iter =10000,n.thin=1,n.burnin = 1000,bugs.directory = "C://Users//atief//Desktop//efe//master//derslerr//2nd semester//bayesian//winbugs//WinBUGS14")
```

``` r
fake_data_model3
```

    ## Inference for Bugs model at "fake_m3.bugs", fit using WinBUGS,
    ##  3 chains, each with 10000 iterations (first 1000 discarded)
    ##  n.sims = 27000 iterations saved
    ##           mean   sd  2.5%   25%   50%   75% 97.5% Rhat n.eff
    ## alpha0     0.0  0.1  -0.2  -0.1   0.0   0.0   0.1    1  3400
    ## alpha1    -0.2  0.1  -0.4  -0.2  -0.2  -0.1   0.1    1 27000
    ## sigma      0.7  0.1   0.6   0.6   0.7   0.7   0.8    1  2400
    ## deviance 901.4 19.8 864.7 887.8 900.8 914.5 942.1    1  8600
    ## 
    ## For each parameter, n.eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
    ## 
    ## DIC info (using the rule, pD = Dbar-Dhat)
    ## pD = 141.2 and DIC = 1042.7
    ## DIC is an estimate of expected predictive error (lower deviance is better).

Again as expected, this much simple model is not as succesfull as the
previos models for describing the data. However it is important to note
that Deviance Information Criteria(DIC) is not a very reliable way to
compare random effect models; hence the results may be misleading.

# Conclusion

In this project, we made a bayesian analysis and inference using Seeds
data. After posterior parameter estimates are obtained, it is compared
to the estimates obtained from frequentist approach. Furthermore, new
data is generated using the estimated parameters. The new generated data
fitted on different models in order to show that it fits best to the
model that created it, as it’s expected intuitively.

# Appendix

Convergence Diagnostics for the models fitted to the new generated data:

Diagnostics about fake\_data\_model1:

``` r
for (i in 1:5){
plot(fake_data_model1$sims.array[,1,i], type = "l", main=names(fake_data_model1$sims.list)[i])
  
}
```

![](project_markdown_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-29-2.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-29-3.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-29-4.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-29-5.png)<!-- -->

``` r
par(mfrow=c(3,2))
for (i in 1:4){
  acf(fake_data_model1$sims.array[,1,i],main=names(fake_data_model1$sims.list)[i])
}
par(mfrow=c(1,1))
```

![](project_markdown_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

Diagnostics about fake\_data\_model2:

``` r
for (i in 1:5){
plot(fake_data_model2$sims.array[,1,i], type = "l", main=names(fake_data_model2$sims.list)[i])
  
}
```

![](project_markdown_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-31-2.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-31-3.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-31-4.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-31-5.png)<!-- -->

``` r
par(mfrow=c(3,2))
for (i in 1:4){
  acf(fake_data_model2$sims.array[,1,i],main=names(fake_data_model2$sims.list)[i])
}
par(mfrow=c(1,1))
```

![](project_markdown_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

Diagnostics for fake\_data\_model3 :

``` r
for (i in 1:4){
plot(fake_data_model3$sims.array[,1,i], type = "l", main=names(fake_data_model3$sims.list)[i])
  
}
```

![](project_markdown_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-33-2.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-33-3.png)<!-- -->![](project_markdown_files/figure-gfm/unnamed-chunk-33-4.png)<!-- -->

``` r
par(mfrow=c(3,2))
for (i in 1:3){
  acf(fake_data_model3$sims.array[,1,i],main=names(fake_data_model3$sims.list)[i])
}
par(mfrow=c(1,1))
```

![](project_markdown_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

Complete Data Generated:

``` r
generated_data_complete=cbind(as.data.frame(n_vals),as.data.frame(k),
                            as.data.frame(x1_vals),as.data.frame(x2_vals)
                            )
names(generated_data_complete)=c('total no', 'germinated', 'x1', 'x2')
print(generated_data_complete)
```

    ##     total no germinated x1 x2
    ## 1         45         20  1  0
    ## 2         42         15  1  1
    ## 3         43         33  0  1
    ## 4         30         21  0  1
    ## 5         14          4  0  0
    ## 6         53         29  1  1
    ## 7         48         37  0  1
    ## 8         13         10  1  1
    ## 9         41         21  1  0
    ## 10        13          5  0  0
    ## 11        31         21  1  1
    ## 12         7          1  1  1
    ## 13        31         20  1  1
    ## 14        15          7  1  0
    ## 15        25         12  1  0
    ## 16        49         13  1  1
    ## 17        44         34  0  1
    ## 18        40         31  0  1
    ## 19        31         13  0  0
    ## 20        28         15  1  1
    ## 21         5          2  1  1
    ## 22         6          2  1  1
    ## 23        40         11  0  0
    ## 24        22          4  0  0
    ## 25        22         11  0  0
    ## 26        25          9  1  0
    ## 27        44         15  0  0
    ## 28         7          5  0  1
    ## 29        50         17  1  0
    ## 30        52         25  1  0
    ## 31         6          6  0  1
    ## 32        20          9  1  1
    ## 33        43         15  1  1
    ## 34        20          5  0  0
    ## 35         6          3  1  0
    ## 36        31         19  0  1
    ## 37         3          2  1  1
    ## 38        24         10  1  0
    ## 39        23          8  0  0
    ## 40        20          8  0  0
    ## 41        19          5  1  0
    ## 42        44         13  1  0
    ## 43        26         17  1  1
    ## 44         8          3  1  1
    ## 45         4          3  1  0
    ## 46        16         12  1  1
    ## 47        15          5  1  0
    ## 48        46         18  1  0
    ## 49        28         16  1  1
    ## 50        43         13  0  0
    ## 51        42         10  0  0
    ## 52        28          8  0  0
    ## 53        45         17  0  0
    ## 54        50         24  1  1
    ## 55        22         10  1  0
    ## 56        37         25  1  1
    ## 57        21          5  0  0
    ## 58        12          2  1  0
    ## 59        42         25  1  1
    ## 60        47         41  0  1
    ## 61        55         36  1  1
    ## 62        38         16  1  0
    ## 63         4          1  0  0
    ## 64        18          6  1  1
    ## 65         5          1  0  0
    ## 66         4          3  0  0
    ## 67         7          3  1  1
    ## 68        53         33  1  1
    ## 69        46         14  0  0
    ## 70        16          4  1  0
    ## 71        38         31  0  1
    ## 72        43         19  0  0
    ## 73        52          7  1  0
    ## 74        32         11  1  0
    ## 75        26         11  1  0
    ## 76        23         16  1  1
    ## 77        49         21  0  0
    ## 78        40         17  0  0
    ## 79        39         15  1  0
    ## 80        38         10  0  0
    ## 81        40         26  1  0
    ## 82        39         14  1  0
    ## 83        22          8  1  1
    ## 84        42         17  0  0
    ## 85        18          7  0  0
    ## 86         7          3  1  0
    ## 87        16         14  0  1
    ## 88        23          8  1  0
    ## 89        46         26  1  1
    ## 90        48         10  1  0
    ## 91        29         15  0  1
    ## 92        18          2  1  0
    ## 93        24         12  1  0
    ## 94        38         32  0  1
    ## 95        46         16  0  0
    ## 96        33         16  1  0
    ## 97        53         24  1  1
    ## 98        12          2  0  0
    ## 99        31         19  1  1
    ## 100       36         12  1  0
    ## 101       21          7  0  0
    ## 102       34          9  0  0
    ## 103       32         23  0  1
    ## 104       12          7  1  0
    ## 105        5          4  0  0
    ## 106       39         19  0  0
    ## 107       45         13  1  0
    ## 108       14          4  0  0
    ## 109       17          6  0  0
    ## 110       28         21  0  1
    ## 111       21         12  1  1
    ## 112       16         12  0  1
    ## 113       21         12  1  1
    ## 114       46         20  1  1
    ## 115       26          7  0  0
    ## 116       28         10  1  0
    ## 117       51         29  1  1
    ## 118        5          3  1  0
    ## 119       19          4  0  0
    ## 120        4          0  0  0
    ## 121        5          1  0  0
    ## 122       27         15  1  1
    ## 123       20         16  0  1
    ## 124       34         13  1  0
    ## 125       46         25  0  1
    ## 126       21          4  0  0
    ## 127       53         12  1  0
    ## 128       23         17  0  1
    ## 129       11          5  1  1
    ## 130       13          1  0  0
    ## 131       53         17  0  0
    ## 132       45         23  0  0
    ## 133       48         29  0  1
    ## 134       20         14  0  1
    ## 135       19         11  1  1
    ## 136       25         12  1  1
    ## 137        9          4  1  1
    ## 138       44         15  0  0
    ## 139       19         11  1  0
    ## 140       51         27  0  1
    ## 141       52         20  1  0
    ## 142       46         37  0  1
    ## 143       13          6  0  1
    ## 144       32         16  0  1
    ## 145       13          3  1  0
    ## 146       45         22  0  0
    ## 147       31         10  0  0
    ## 148       50         18  0  0
    ## 149        6          1  1  0
    ## 150       55         27  1  0
    ## 151       51         24  0  1
    ## 152       19         15  1  1
    ## 153       52         18  1  0
    ## 154       39         22  0  1
    ## 155       11          1  1  0
    ## 156        3          1  0  0
    ## 157       49         33  0  1
    ## 158       19          5  0  0
    ## 159       10          4  0  0
    ## 160       25          5  1  0
    ## 161       12          7  1  0
    ## 162       49         35  0  1
    ## 163       44         20  0  1
    ## 164       50         32  1  1
    ## 165       20          6  1  0
    ## 166       42         15  1  0
    ## 167        5          4  1  1
    ## 168        7          3  1  0
    ## 169       32         22  0  1
    ## 170       49         28  0  1
    ## 171       45         23  1  1
    ## 172       14          6  1  0
    ## 173       24         12  0  0
    ## 174       30         19  1  1
    ## 175       28         11  1  0
    ## 176        9          2  1  0
    ## 177       41         25  1  1
    ## 178       24          6  1  0
    ## 179       17          2  0  0
    ## 180        8          3  1  1
    ## 181       53         35  1  1
    ## 182       39         29  0  1
    ## 183       17         10  1  1
    ## 184       35         19  0  1
    ## 185       46         37  1  1
    ## 186       47          9  1  0
    ## 187       37         25  0  1
    ## 188       31         12  0  0
    ## 189       55         34  1  1
    ## 190       26          9  1  0
    ## 191        4          4  0  1
    ## 192       38         14  1  0
    ## 193       27         18  0  0
    ## 194        9          4  0  0
    ## 195       17          7  0  0
    ## 196       45         27  0  1
    ## 197       20         12  1  0
    ## 198       53         16  1  1
    ## 199       28         11  0  0
    ## 200       51         24  1  1
