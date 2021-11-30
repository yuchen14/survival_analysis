Survival Analysis of Adjunct Therapy for Colon Cancer
================
Yuchen Zheng
11/24/2020

``` r
library(survminer)
```

    ## Loading required package: ggplot2

    ## Loading required package: ggpubr

``` r
library(ggplot2)
library(survival)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
#convert numberic variables to factors 
colon$sex <- factor(colon$sex) 
colon$obstruct <- as.factor(colon$obstruct) 
colon$perfor<- factor(colon$perfor) 
colon$adhere <- factor(colon$adhere) 
colon$differ <- factor(colon$differ)
colon$extent <- factor(colon$extent) 
colon$surg <- factor(colon$surg) 
colon$node4 <- factor(colon$node4) 
colon$etype<- factor(colon$etype)
```

``` r
#subset recurrence data
colon.recurrence <- subset(colon, etype == 1, select=c(id:etype))
```

``` r
#plot KM estimate for recurrence data
recurrence.surv <- Surv(colon.recurrence$time, colon.recurrence$status)
recurrence.fit <- survfit(recurrence.surv~rx, data=colon.recurrence)
ggsurvplot(recurrence.fit, conf.int=F,
           title = 'Kaplan-Meier Curve for Colon Cancer \nRecurrence by Treatment',
           xlab= 'Times (until death) \n in Days')
```

![](Survival-Analysis_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
r.rx.coxph <- coxph(Surv(time, status) ~ rx, data=colon.recurrence)
summary(r.rx.coxph)
```

    ## Call:
    ## coxph(formula = Surv(time, status) ~ rx, data = colon.recurrence)
    ## 
    ##   n= 929, number of events= 468 
    ## 
    ##               coef exp(coef) se(coef)      z Pr(>|z|)    
    ## rxLev     -0.01512   0.98499  0.10708 -0.141    0.888    
    ## rxLev+5FU -0.51209   0.59924  0.11863 -4.317 1.58e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##           exp(coef) exp(-coef) lower .95 upper .95
    ## rxLev        0.9850      1.015    0.7985    1.2150
    ## rxLev+5FU    0.5992      1.669    0.4749    0.7561
    ## 
    ## Concordance= 0.554  (se = 0.013 )
    ## Likelihood ratio test= 24.34  on 2 df,   p=5e-06
    ## Wald test            = 22.58  on 2 df,   p=1e-05
    ## Score (logrank) test = 23.07  on 2 df,   p=1e-05

``` r
#removing NA values
colon.recurrence1 = na.omit(colon.recurrence)
```

``` r
r.model1 <- coxph(Surv(time, status) ~ sex, data=colon.recurrence1)
r.model2 <- coxph(Surv(time, status) ~ age, data=colon.recurrence1)
r.model3 <- coxph(Surv(time, status) ~ obstruct, data=colon.recurrence1)
r.model4 <- coxph(Surv(time, status) ~ perfor, data=colon.recurrence1)
r.model5 <- coxph(Surv(time, status) ~ adhere, data=colon.recurrence1)
r.model6 <- coxph(Surv(time, status) ~ differ, data=colon.recurrence1)
r.model7 <- coxph(Surv(time, status) ~ extent, data=colon.recurrence1)
r.model8 <- coxph(Surv(time, status) ~ surg, data=colon.recurrence1)
r.model9 <- coxph(Surv(time, status) ~ node4, data=colon.recurrence1)
BIC(r.model1, r.model2,r.model3, r.model4,r.model5,r.model6,r.model7,r.model8,r.model9)
```

    ##          df      BIC
    ## r.model1  1 5758.946
    ## r.model2  1 5758.128
    ## r.model3  1 5756.859
    ## r.model4  1 5758.395
    ## r.model5  1 5755.234
    ## r.model6  2 5756.311
    ## r.model7  3 5744.420
    ## r.model8  1 5755.739
    ## r.model9  1 5689.597

``` r
r.model9.1 <- coxph(Surv(time, status) ~ node4+sex, data=colon.recurrence1)
r.model9.2 <- coxph(Surv(time, status) ~ node4+age, data=colon.recurrence1)
r.model9.3 <- coxph(Surv(time, status) ~ node4+obstruct, data=colon.recurrence1)
r.model9.4 <- coxph(Surv(time, status) ~ node4+perfor, data=colon.recurrence1)
r.model9.5 <- coxph(Surv(time, status) ~ node4+adhere, data=colon.recurrence1)
r.model9.6 <- coxph(Surv(time, status) ~ node4+differ, data=colon.recurrence1)
r.model9.7 <- coxph(Surv(time, status) ~ node4+extent, data=colon.recurrence1)
r.model9.8 <- coxph(Surv(time, status) ~ node4+surg, data=colon.recurrence1)
BIC(r.model9.1,r.model9.2,r.model9.3,r.model9.4,r.model9.5,r.model9.6,r.model9.7,r.model9.8)
```

    ##            df      BIC
    ## r.model9.1  2 5694.438
    ## r.model9.2  2 5695.111
    ## r.model9.3  2 5690.740
    ## r.model9.4  2 5693.106
    ## r.model9.5  2 5691.062
    ## r.model9.6  3 5695.820
    ## r.model9.7  4 5685.235
    ## r.model9.8  2 5689.928

``` r
r.model9.7.1 <- coxph(Surv(time, status) ~ node4+extent+sex, data=colon.recurrence1)
r.model9.7.2 <- coxph(Surv(time, status) ~ node4+extent+age, data=colon.recurrence1)
r.model9.7.3 <- coxph(Surv(time, status) ~ node4+extent+obstruct, data=colon.recurrence1)
r.model9.7.4 <- coxph(Surv(time, status) ~ node4+extent+perfor, data=colon.recurrence1)
r.model9.7.5 <- coxph(Surv(time, status) ~ node4+extent+adhere, data=colon.recurrence1)
r.model9.7.6 <- coxph(Surv(time, status) ~ node4+extent+differ, data=colon.recurrence1)
r.model9.7.7 <- coxph(Surv(time, status) ~ node4+extent+surg, data=colon.recurrence1)
BIC(r.model9.7.1,r.model9.7.2,r.model9.7.3, r.model9.7.4,r.model9.7.5,r.model9.7.6,r.model9.7.7)
```

    ##              df      BIC
    ## r.model9.7.1  5 5690.331
    ## r.model9.7.2  5 5690.806
    ## r.model9.7.3  5 5688.145
    ## r.model9.7.4  5 5689.742
    ## r.model9.7.5  5 5688.623
    ## r.model9.7.6  6 5692.914
    ## r.model9.7.7  5 5685.198

``` r
r.model9.7.7.1 <- coxph(Surv(time, status) ~ node4+extent+surg+sex, data=colon.recurrence1)
r.model9.7.7.2 <- coxph(Surv(time, status) ~ node4+extent+surg+age, data=colon.recurrence1)
r.model9.7.7.3 <- coxph(Surv(time, status) ~ node4+extent+surg+obstruct, data=colon.recurrence1)
r.model9.7.7.4 <- coxph(Surv(time, status) ~ node4+extent+surg+perfor, data=colon.recurrence1)
r.model9.7.7.5 <- coxph(Surv(time, status) ~ node4+extent+surg+adhere, data=colon.recurrence1)
r.model9.7.7.6 <- coxph(Surv(time, status) ~ node4+extent+surg+differ, data=colon.recurrence1)
BIC(r.model9.7.7.1, r.model9.7.7.2,r.model9.7.7.3,r.model9.7.7.4,r.model9.7.7.5,r.model9.7.7.6)
```

    ##                df      BIC
    ## r.model9.7.7.1  6 5690.164
    ## r.model9.7.7.2  6 5690.688
    ## r.model9.7.7.3  6 5688.372
    ## r.model9.7.7.4  6 5689.722
    ## r.model9.7.7.5  6 5688.796
    ## r.model9.7.7.6  7 5693.056

``` r
#full model
r.model.full <- coxph(Surv(time, status) ~ sex + age + obstruct + perfor + adhere + differ + extent + surg + node4, data=colon.recurrence1)
```

``` r
BIC(r.model.full,r.model9.7.7.3,r.model9.7.7,r.model9.7,r.model9)
```

    ##                df      BIC
    ## r.model.full   12 5716.259
    ## r.model9.7.7.3  6 5688.372
    ## r.model9.7.7    5 5685.198
    ## r.model9.7      4 5685.235
    ## r.model9        1 5689.597

We can use the Analysis of Deviance procedure to get the proper
Likelihood Ratio Test

``` r
anova(r.model9.7.7)
```

    ## Analysis of Deviance Table
    ##  Cox model: response is Surv(time, status)
    ## Terms added sequentially (first to last)
    ## 
    ##         loglik   Chisq Df Pr(>|Chi|)    
    ## NULL   -2877.1                          
    ## node4  -2841.8 70.7685  1  < 2.2e-16 ***
    ## extent -2830.4 22.6631  3  4.747e-05 ***
    ## surg   -2827.3  6.1378  1    0.01323 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
r.model9.7.7
```

    ## Call:
    ## coxph(formula = Surv(time, status) ~ node4 + extent + surg, data = colon.recurrence1)
    ## 
    ##            coef exp(coef) se(coef)     z       p
    ## node41  0.84754   2.33389  0.09905 8.557 < 2e-16
    ## extent2 0.31142   1.36537  0.52970 0.588 0.55658
    ## extent3 0.88437   2.42147  0.50427 1.754 0.07947
    ## extent4 1.44544   4.24373  0.54141 2.670 0.00759
    ## surg1   0.26138   1.29872  0.10361 2.523 0.01165
    ## 
    ## Likelihood ratio test=99.57  on 5 df, p=< 2.2e-16
    ## n= 888, number of events= 446

``` r
colon.recurrence2 <- colon.recurrence
colon.recurrence2$nodes[is.na(colon.recurrence2$nodes)] <- mean(colon.recurrence2$nodes, na.rm = TRUE)
colon.recurrence2$differ[is.na(colon.recurrence2$differ)] <- 2
colon.recurrence2$differ <- factor(colon.recurrence2$differ, exclude=NULL)
```

We created to a second recurrence dataset by replcing the na values in
nodes with the mean of the column and include the NA value as a
character level in the differ column to check the cox proportional model
assumption.

``` r
recurrence2.surv <- Surv(colon.recurrence2$time, colon.recurrence2$status)
plot(survfit(recurrence2.surv~rx, data=colon.recurrence2),
     fun='cloglog',
     col =c(2,3,4),
     ylab='log-log(S(t))',
     xlab='Time(days)',
     main="C-log-log plot for Covariate rx")
legend("topleft" ,legend =c("Obs","Lev","Lev+5FU"), col=c(2,3,4),lty=1)
```

![](Survival-Analysis_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->
The three curves in the cloglog plot for covariate rx is overlapping
before 100 days, but after 100 days the two curves for Obs and Lev are
crossing over and the two curves appear to be parallel to the curve for
Lev+5Fu. We are somewhat concerning that the assumption might be
violated because of the crossover before 100 days.

``` r
plot(survfit(recurrence2.surv~surg, data=colon.recurrence2),
     fun='cloglog',
     col =c(2,3),
     ylab='log-log(S(t))',
     xlab='Time(days)',
     main="C-log-log plot for Covariate surg")
legend("topleft" ,legend =c("short","long"), col=c(2,3),lty=1)
```

![](Survival-Analysis_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->
The two curves in the C-log-log plot for covariate surg is crossing over
and very close to each other. We conclude that the assumption is
approriate to use for this covariate because crossingover means the two
groups have the same hazard ratio.

``` r
plot(survfit(recurrence2.surv~extent, data=colon.recurrence2),
     fun='cloglog',
     col =c(2,3,4,5),
     ylab='log-log(S(t))',
     xlab='Time(days)',
     main="C-log-log plot for Covariate extent")
legend("topleft" ,legend =c("1-submucosa","2-muscle","3-serosa", "4-contiguous structures"), col=c(2,3,4,5),lty=1)
```

![](Survival-Analysis_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->
The curves in the C-log-log plot are crossing over after 100 days. Since
there are not enough data points in each group to show a more
comprehesive trend, it’s hard for us to make a decision based on the
plot.

``` r
plot(survfit(recurrence2.surv~node4, data=colon.recurrence2),
     fun='cloglog',
     col =c(2,3),
     ylab='log-log(S(t))',
     xlab='Time(days)',
     main="C-log-log plot for Covariate node4")
legend("topleft" ,legend =c("more than 4 postivie nodes","less than 4 nodes"), col=c(2,3),lty=1)
```

![](Survival-Analysis_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->
The two curves in this C-log-log plot appear to be parallel to each
other, so we think the cox proportional assumption is approritate for
the covariate
node4.

``` r
cox.zph(coxph(formula = Surv(time, status) ~ node4 +extent+surg+rx, data = colon.recurrence2))
```

    ##          chisq df       p
    ## node4  10.9686  1 0.00093
    ## extent  1.6141  3 0.65620
    ## surg    1.2254  1 0.26830
    ## rx      0.0866  2 0.95761
    ## GLOBAL 14.1390  7 0.04876

Since the p value for node4 is less than 0.05, there is significant
evidence that the cox proportional model assumption is violated for
variable node4. However, the C-log-log plot for covariate node4 looks
parallel, so we think the significance of p-value might be due to some
noise in the data at the beginning of the study. Additionally, the
overall p-value is greater than 0.05 which means we fail to reject the
null hypothesis and can conclude that the cox proportional model
assumption is reasonable to use for this
model.

``` r
colon.recurrence.coxph <- coxph(Surv(time, status) ~ node4+extent++surg+rx, data=colon.recurrence1)
summary(colon.recurrence.coxph)
```

    ## Call:
    ## coxph(formula = Surv(time, status) ~ node4 + extent + +surg + 
    ##     rx, data = colon.recurrence1)
    ## 
    ##   n= 888, number of events= 446 
    ## 
    ##               coef exp(coef) se(coef)      z Pr(>|z|)    
    ## node41     0.84049   2.31750  0.09908  8.483  < 2e-16 ***
    ## extent2    0.26290   1.30070  0.53001  0.496   0.6199    
    ## extent3    0.84593   2.33014  0.50468  1.676   0.0937 .  
    ## extent4    1.38116   3.97953  0.54157  2.550   0.0108 *  
    ## surg1      0.23556   1.26562  0.10399  2.265   0.0235 *  
    ## rxLev     -0.02871   0.97170  0.11021 -0.260   0.7945    
    ## rxLev+5FU -0.49312   0.61072  0.12169 -4.052 5.08e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##           exp(coef) exp(-coef) lower .95 upper .95
    ## node41       2.3175     0.4315    1.9085    2.8142
    ## extent2      1.3007     0.7688    0.4603    3.6756
    ## extent3      2.3301     0.4292    0.8666    6.2656
    ## extent4      3.9795     0.2513    1.3767   11.5032
    ## surg1        1.2656     0.7901    1.0323    1.5517
    ## rxLev        0.9717     1.0291    0.7829    1.2060
    ## rxLev+5FU    0.6107     1.6374    0.4811    0.7752
    ## 
    ## Concordance= 0.656  (se = 0.013 )
    ## Likelihood ratio test= 120.3  on 7 df,   p=<2e-16
    ## Wald test            = 122.2  on 7 df,   p=<2e-16
    ## Score (logrank) test = 129.6  on 7 df,   p=<2e-16

``` r
anova(colon.recurrence.coxph)
```

    ## Analysis of Deviance Table
    ##  Cox model: response is Surv(time, status)
    ## Terms added sequentially (first to last)
    ## 
    ##         loglik   Chisq Df Pr(>|Chi|)    
    ## NULL   -2877.1                          
    ## node4  -2841.8 70.7685  1  < 2.2e-16 ***
    ## extent -2830.4 22.6631  3  4.747e-05 ***
    ## surg   -2827.3  6.1378  1    0.01323 *  
    ## rx     -2817.0 20.7443  2  3.129e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

\#Marginal model for Death

``` r
colon.death <- subset(colon, etype == 2)

death.fit <- survfit(Surv(time,status) ~ rx, data = colon.death)
ggsurvplot(death.fit, conf.int = F, 
           title = "Kaplan-Meier Curve for Colon Cancer Mortality \nby Treatment", 
           xlab = "Time (until death) \n in Days")
```

![](Survival-Analysis_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
d.rx.coxph <- coxph(Surv(time, status) ~ rx, data=colon.death)
summary(d.rx.coxph)
```

    ## Call:
    ## coxph(formula = Surv(time, status) ~ rx, data = colon.death)
    ## 
    ##   n= 929, number of events= 452 
    ## 
    ##               coef exp(coef) se(coef)      z Pr(>|z|)   
    ## rxLev     -0.02664   0.97371  0.11030 -0.241  0.80917   
    ## rxLev+5FU -0.37171   0.68955  0.11875 -3.130  0.00175 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##           exp(coef) exp(-coef) lower .95 upper .95
    ## rxLev        0.9737      1.027    0.7844    1.2087
    ## rxLev+5FU    0.6896      1.450    0.5464    0.8703
    ## 
    ## Concordance= 0.536  (se = 0.013 )
    ## Likelihood ratio test= 12.15  on 2 df,   p=0.002
    ## Wald test            = 11.56  on 2 df,   p=0.003
    ## Score (logrank) test = 11.68  on 2 df,   p=0.003

``` r
#removing NA values
colon.death1 <- na.omit(colon.death)
```

``` r
d.model1 <- coxph(Surv(time, status) ~ sex, data=colon.death1)
d.model2 <- coxph(Surv(time, status) ~ age, data=colon.death1)
d.model3 <- coxph(Surv(time, status) ~ obstruct, data=colon.death1)
d.model4 <- coxph(Surv(time, status) ~ perfor, data=colon.death1)
d.model5 <- coxph(Surv(time, status) ~ adhere, data=colon.death1)
d.model6 <- coxph(Surv(time, status) ~ node4, data=colon.death1)
d.model7 <- coxph(Surv(time, status) ~ differ, data=colon.death1)
d.model8 <- coxph(Surv(time, status) ~ extent, data=colon.death1)
d.model9 <- coxph(Surv(time, status) ~ surg, data=colon.death1)


BIC(d.model1, d.model2, d.model3, d.model4, d.model5, d.model6, d.model7, d.model8, d.model9)
```

    ##          df      BIC
    ## d.model1  1 5541.880
    ## d.model2  1 5541.408
    ## d.model3  1 5537.618
    ## d.model4  1 5541.525
    ## d.model5  1 5536.616
    ## d.model6  1 5458.201
    ## d.model7  2 5535.900
    ## d.model8  3 5526.573
    ## d.model9  1 5538.031

``` r
d.model6.1 <- coxph(Surv(time, status) ~ node4 + sex, data=colon.death1)
d.model6.2 <- coxph(Surv(time, status) ~ node4 + age, data=colon.death1)
d.model6.3 <- coxph(Surv(time, status) ~ node4 + obstruct, data=colon.death1)
d.model6.4 <- coxph(Surv(time, status) ~ node4 + perfor, data=colon.death1)
d.model6.5 <- coxph(Surv(time, status) ~ node4 + adhere, data=colon.death1)
d.model6.6 <- coxph(Surv(time, status) ~ node4 + surg, data=colon.death1)
d.model6.7 <- coxph(Surv(time, status) ~ node4 + differ, data=colon.death1)
d.model6.8 <- coxph(Surv(time, status) ~ node4 + extent, data=colon.death1)

BIC(d.model6.1, d.model6.2, d.model6.3, d.model6.4, d.model6.5, d.model6.6, d.model6.7, d.model6.8)
```

    ##            df      BIC
    ## d.model6.1  2 5464.237
    ## d.model6.2  2 5461.289
    ## d.model6.3  2 5458.564
    ## d.model6.4  2 5463.690
    ## d.model6.5  2 5459.417
    ## d.model6.6  2 5458.843
    ## d.model6.7  3 5463.982
    ## d.model6.8  4 5455.598

``` r
d.model6.8.1 <- coxph(Surv(time, status) ~ node4 + extent + sex, data=colon.death1)
d.model6.8.2 <- coxph(Surv(time, status) ~ node4 + extent + age, data=colon.death1)
d.model6.8.3 <- coxph(Surv(time, status) ~ node4 + extent + obstruct, data=colon.death1)
d.model6.8.4 <- coxph(Surv(time, status) ~ node4 + extent + perfor, data=colon.death1)
d.model6.8.5 <- coxph(Surv(time, status) ~ node4 + extent + adhere, data=colon.death1)
d.model6.8.6 <- coxph(Surv(time, status) ~ node4 + extent + differ, data=colon.death1)
d.model6.8.7 <- coxph(Surv(time, status) ~ node4 + extent + surg, data=colon.death1)

BIC(d.model6.8.1, d.model6.8.2, d.model6.8.3, d.model6.8.4, d.model6.8.5, d.model6.8.6, d.model6.8.7)
```

    ##              df      BIC
    ## d.model6.8.1  5 5461.607
    ## d.model6.8.2  5 5458.739
    ## d.model6.8.3  5 5457.686
    ## d.model6.8.4  5 5461.424
    ## d.model6.8.5  5 5458.507
    ## d.model6.8.6  6 5462.390
    ## d.model6.8.7  5 5456.028

``` r
d.model.full <- coxph(Surv(time, status) ~ sex + age + obstruct + perfor + adhere + differ + extent + surg + node4, data=colon.death1)
```

``` r
BIC(d.model.full,d.model6.8.7, d.model6.8, d.model6)
```

    ##              df      BIC
    ## d.model.full 12 5483.868
    ## d.model6.8.7  5 5456.028
    ## d.model6.8    4 5455.598
    ## d.model6      1 5458.201

``` r
anova(d.model6.8)
```

    ## Analysis of Deviance Table
    ##  Cox model: response is Surv(time, status)
    ## Terms added sequentially (first to last)
    ## 
    ##         loglik  Chisq Df Pr(>|Chi|)    
    ## NULL   -2767.9                         
    ## node4  -2726.1 83.685  1  < 2.2e-16 ***
    ## extent -2715.7 20.794  3  0.0001162 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
colon.death2 <- colon.death
colon.death2$nodes[is.na(colon.death2$nodes)] <- mean(colon.death2$nodes, na.rm = TRUE)

colon.death2$differ <- factor(colon.death2$differ, exclude=NULL)
```

``` r
dnode4.fit <- survfit(Surv(time, status) ~ node4, data = colon.death2)

ggsurvplot(dnode4.fit , conf.int = F, 
           fun = "cloglog",
           xlim = c(20, 5000),
           title = "C-Log-Log for Colon Cancer Mortality \nby node4", 
           xlab = "Time (until death) \n in Days")
```

![](Survival-Analysis_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

``` r
dextent.fit <- survfit(Surv(time, status) ~ extent, data = colon.death2)

ggsurvplot(dextent.fit , conf.int = F, 
           fun = "cloglog",
           xlim = c(20, 5000),
           title = "C-Log-Log for Colon Cancer Mortality \nby extent", 
           xlab = "Time (until death) \n in Days")
```

![](Survival-Analysis_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

``` r
dextent.fit <- survfit(Surv(time, status) ~ rx, data = colon.death2)

ggsurvplot(dextent.fit , conf.int = F, 
           fun = "cloglog",
           xlim = c(20, 5000),
           title = "C-Log-Log for Colon Cancer Mortality \nby rx", 
           xlab = "Time (until death) \n in Days")
```

![](Survival-Analysis_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

``` r
cox.zph(coxph(formula = Surv(time, status) ~ node4 +extent+rx, data = colon.death2))
```

    ##        chisq df     p
    ## node4   6.50  1 0.011
    ## extent  7.24  3 0.065
    ## rx      1.99  2 0.369
    ## GLOBAL 14.79  6 0.022

``` r
d.final.model <- coxph(Surv(time, status) ~ node4 + extent + rx, data=colon.death2)
summary(d.final.model)
```

    ## Call:
    ## coxph(formula = Surv(time, status) ~ node4 + extent + rx, data = colon.death2)
    ## 
    ##   n= 929, number of events= 452 
    ## 
    ##               coef exp(coef) se(coef)      z Pr(>|z|)    
    ## node41     0.91540   2.49779  0.09702  9.435  < 2e-16 ***
    ## extent2    0.41635   1.51642  0.52816  0.788  0.43051    
    ## extent3    0.93099   2.53703  0.50453  1.845  0.06500 .  
    ## extent4    1.41270   4.10702  0.53499  2.641  0.00828 ** 
    ## rxLev     -0.04248   0.95841  0.11052 -0.384  0.70071    
    ## rxLev+5FU -0.37917   0.68443  0.11890 -3.189  0.00143 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##           exp(coef) exp(-coef) lower .95 upper .95
    ## node41       2.4978     0.4004    2.0652     3.021
    ## extent2      1.5164     0.6594    0.5386     4.270
    ## extent3      2.5370     0.3942    0.9438     6.820
    ## extent4      4.1070     0.2435    1.4393    11.720
    ## rxLev        0.9584     1.0434    0.7718     1.190
    ## rxLev+5FU    0.6844     1.4611    0.5422     0.864
    ## 
    ## Concordance= 0.651  (se = 0.013 )
    ## Likelihood ratio test= 122.5  on 6 df,   p=<2e-16
    ## Wald test            = 126.1  on 6 df,   p=<2e-16
    ## Score (logrank) test = 135.8  on 6 df,   p=<2e-16

``` r
anova(d.final.model)
```

    ## Analysis of Deviance Table
    ##  Cox model: response is Surv(time, status)
    ## Terms added sequentially (first to last)
    ## 
    ##         loglik  Chisq Df Pr(>|Chi|)    
    ## NULL   -2930.2                         
    ## node4  -2885.7 88.954  1  < 2.2e-16 ***
    ## extent -2875.1 21.321  3   9.03e-05 ***
    ## rx     -2869.0 12.177  2   0.002269 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
