# ----------------------------------------------------------------- Introduction -----------------------------------------------------------
#Many problems in both science and engineering require fitting a distributional model to a
#uni-variate dataset. That is, data consisting of a set of N empirical observations obtained
#by measuring a certain system, provided that the measurements are independent and come
#from the same distribution. Distributional modeling can be used in several contexts. Statistical tests typically depend 
#on certain assumptions with regards to the underlying distribution (e.g., many tests are based on the assumption of normality). 
#Appropriate distributional models are also usedto more accurately assess uncertainty and, in particular, to assess the uncertainty over the
#full range of the data. Although a poorly chosen distributional model may suffice for measuring and assessing the uncertainty of averages, 
#this will not be the case for tails of the distribution. In many applications (e.g., reliability), accurate assessment of the tail behavior
#is more critical than the average. Obtaining an appropriate distributional model can be an
#exhaustive process that takes time, patience and requires previous knowledge of statistics
#and is, therefore, a difficult task for some analysts.






# -------------------------------------------------- Steps included in finding a statistical distribution --------------------------------------
# 1. Facilitate the assumptions of tests to be done ()
# 1. Create/find some data - hopefully many
# 2. Vizualise the sample/kernel distribution to get undestanding of what type/family/continous/discrete its generated from.
# 3. Use some sort of goodness-of-fit test to compare this sample data to one or more sample distributions.
# 4. Repeat step 3 until satisfactory.

# ------------------------------------------------------------- Creation of "random sample" ----------------------------------------------------

SampleNorm34 <- rexp(10000,11)                                              # Lets create an "unknown" sample of X ~ N(3,4)
Dnorm <- function(x,my,sd){                                                 # Also create a function with manually written normal dist formula 
 1 / (sqrt(2*pi)*sd) * exp(-((x-my)^2)/(2*sd^2))                            # to show that it is equal to plot it as with "dnorm" but also to show
}                                                                           # that plot(dnorm(x,3,4)) does not work.


# ---------------------------------------------------------------- Follow Assumptions -----------------------------------------------------------
# 1. The null hypothesis is both samples are randomly drawn from the same (pooled) set of values.
# 2. The two samples are mutually independent.
# 3. The scale of measurement is at least ordinal.
# 4. The test is only exact for continuous variables. It is conservative for discrete variables.

# First of we see that the values in SampleNorm34 are both at least ordinal and continuous. This first point is what we will check later BUT
# we can check for independence. We will use chi-square test for independence.
chisq.test(SampleNorm34, rnorm(10000,3,4))
# Since this is the closes we come to our sample dist AND STILL INDEPENDENT we dont have to do it for other distributions or parameter values.
# By using rnorm etc we always get randomly independent distributed values.



# ----------------------------------------------------- First Visualize the sample distribution --------------------------------------------------
# By looking at the plots we see that its normally distributed with a certain mean. Lets however ignore this and say that we dont know
# anything and that we need to look for all continuous distributions for many parameter values.

plot(SampleNorm34, dnorm(SampleNorm34,3,4))                                # Produces same plot of dist, just to make the point.
plot(SampleNorm34, Dnorm(SampleNorm34,3,4))
plot(density(SampleNorm34))                                                # Plot the kernel denisty (sample density) can be good to compare with plot(x)


# ----------------------------------------------------- First Vizualise the sample distribution ---------------------------------------------------
# Two-sample KS test since we also need to compute correct parameter values (not only distribution). Using a one-sample test would only check
# if our sample belong to that theoretical distribution. By including a large random sample with 100 combinations of discrete parameter values we may
# also find what parameter values might suit our sample distribution. This is why we have multiple candidate distributions below. One could say that
# this looping skips the MLE/MoM to find likely parameter estimates. Below are parameter estimate loops which goes
# 1,1 1,2 1,3 ... 99,97 99,98 99,99 100,1 100,2 ... 100,98 100,99 100,100 and all major continuous statistical distributions. NOTE that the p-value
# in the "IF"-statement will not exactly be the p-value saved in "vec" because of everytime rnorm etc is run it generates new random value. To mitigate
# this we have a set.seed() in each loop to replicate the random values from rnorm etc to generate the same p-value.



GOFTest <- function(X,alpha,n,nparam,prob='No'){

# Create vectors of all n parameter discrete/continuous dist to loop through and try all combos for all dists.  
  DiscDist0ProbParam <- c(rgeom)
  DiscDist1ProbParam <- c(rbinom)
  DiscDist1Param <- c(rpois)
  ContDist0Param <- c(rchisq,rt)
  ContDist2Param <- c(rbeta,rcauchy,rgamma,rlogis,rlnorm,rnorm,rweibull)   # No need for rexp since special case of gamma(1,n)

# Here we state all dists that are continous and have 2 params
  if (typeof(X) == 'double' & nparam == 2){
    k <- 1                                                                 # Create a scalar k to keep track of loop and to what dist in dist vector gives any result
    for (param in ContDist2Param) {
      vec <- c()                                                           # Here, ks.test()$p.value gives a p-value, so that vec will always contain a p-value. The vec will reset each time a new top loop p begins
      for (i in 1:50) {
        for (j in 1:50) {
          set.seed(123)
          if (ks.test(x = X, y = param(n,i,j), alternative = 'two.sided')$p.value > 2*alpha){   # Must have saved vec before printing vec and iterations since rnorm will create new print p-value otherwise. 2*alpha = two-sided CI
            vec <- c(vec, ks.test(x = X, y = param(n,i,j), alternative = 'two.sided')$p.value)  # Store the previous vec value in new vec which .append(x) the values.
            print(paste("p-value =",round(vec,3), "par1 =",i, "par2 =",j, "dist=",ContDist2Param[k]))  # Print p-value, parameters and dist[k] which gives kth element in the vector - which is the correct dist for the data
          } 
         }
        }
       k <- k+1                                                            # For each top loop k = k+1 and saves new k value to use when displaying correct dist for data.
      }
    
# Here we state all dists that are continous and have 0 params (only n & df=n)
  } else if (typeof(X) == 'double' & nparam == 0) {
    k <- 1                                                                 # Re-start counter k outside loop
    for (param in ContDist0Param) {
      vec <- c()
      for (N in 1:50) {
        set.seed(123)
        if (ks.test(x = X, y = param(n,N), alternative = 'two.sided')$p.value > 2*alpha){ 
          vec <- c(vec, ks.test(x = X, y = param(n,N), alternative = 'two.sided')$p.value)
          print(paste("p-value =",round(vec,3), "par1 =",N, "dist=",ContDist0Param[k]))
        }
       }
      k <- k+1
     }
  
# Here we state all dists that are discrete and have 1 params
 } else if (typeof(X) == 'integer' & nparam == 1) {
  k <- 1                                                                   # Re-start counter k outside loop
  for (param in DiscDist1Param) {
    vec <- c()
    for (i in 1:50) {
      set.seed(123)
      if (ks.test(x = X, y = param(n,i), alternative = 'two.sided')$p.value > 2*alpha){ 
        vec <- c(vec, ks.test(x = X, y = param(n,i), alternative = 'two.sided')$p.value)
        print(paste("p-value =",round(vec,3), "par1 =",i, "dist=",DiscDist1Param[k]))
      }
     }
    k <- k+1
   }
  
  # Here we state all dists that are discrete and have 1 param and 1 prob
 } else if (typeof(X) == 'integer' & nparam == 1 & prob == 'Yes'){
    k <- 1                                                                      
    for (param in DiscDist1ProbParam) {
      vec <- c()                                                             
      for (i in 1:50) {
        for (j in seq(0.001,1,len=50)) {                                   # Cant start loop with 0 prob
          set.seed(123)
          if (ks.test(x = X, y = param(n,i,j), alternative = 'two.sided')$p.value > 2*alpha){   
            vec <- c(vec, ks.test(x = X, y = param(n,i,j), alternative = 'two.sided')$p.value)  
            print(paste("p-value =",round(vec,3), "par1 =",i, "par2 =",j, "dist=",DiscDist1ProbParam[k]))  
          } 
         }
        }
       k <- k+1                                                                
      }
  
# Here we state all dists that are discrete and have 0 params (prob)
 } else if (typeof(X) == 'integer' & nparam == 0 & prob == 'Yes') {
  k <- 1                                                                   # Re-start counter k outside loop
  for (param in DiscDist0ProbParam) {
    vec <- c()
    for (i in seq(0.001,1,len=50)) {                                       # Cant start loop with 0 prob
      set.seed(123)
      if (ks.test(x = X, y = param(n,i), alternative = 'two.sided')$p.value > 2*alpha){ 
        vec <- c(vec, ks.test(x = X, y = param(n,i), alternative = 'two.sided')$p.value)
        print(paste("p-value =",round(vec,3), "par1 =",i, "dist=",DiscDist0ProbParam[k]))
      }
     }
    k <- k+1
   }
 } else {print("Please follow input constraints")}
}




# Since no Goodness-of-fit test returns any parameter values except the normal distribution tells us that our sample data was 
# generated from this distribution. The iteration i and j tells us what parameter values was most likely and minimized residuals
# for our theoretical distribution. 
# ANSWER: Xi ~ N(3,4)
# PS, if we take choose continuous values 0-3 for my and 0-4 for sigma, we only get significant GOF-values starting from 2.95 for my and 3.95 for sigma
# further strengthening our argument for N(3,4).

# PS 2: to justify that our sample data has parameters mean:3 and sd:4 we use MLE

library(stats4)

fun <- function(my, sigma){                                                # State function to optimize and its "unknown parameters" - dont use xi=SampleNorm34 since its global and known
  samp <- dnorm(x=SampleNorm34, mean=my, sd=abs(sigma), log = TRUE)        # MUST know underlying distribution BEFOREHAND and must use our sample as xi. MUST have abs(sigma) since no negative variances exist and algorithm doesnt care    
  return(-sum(samp))                                                       
}
(MLE <- mle(minuslogl=fun, start=list(my=0,sigma=1), method="L-BFGS-B"))        # Here we go my ~ 2.979 and sigma ~ 3.996          

# Please try the GOFtest function for perhaps a Weibull dist or something fun.
