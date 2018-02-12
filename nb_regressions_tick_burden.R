
library(MASS)
LSUM<-t2$LSUM
Sex<-as.numeric(t2$Sex)
Treatment<-t2$Treatment
kruskal.test(LSUM ~ Treatment)
NB0<-glm.nb(LSUM ~ Treatment)
NB1<-glm.nb(LSUM ~ Sex + Treatment)
NB2<-glm.nb(LSUM ~ Sex + Treatment + Sex*Treatment)

 summary(NB0)

Call:
glm.nb(formula = LSUM ~ Treatment, init.theta = 1.388711224, 
    link = log)

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-2.3186  -1.0143  -0.3248   0.2202   2.8177  

Coefficients:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)  2.10802    0.06586  32.006   <2e-16 ***
TreatmentX  -0.22575    0.09449  -2.389   0.0169 *  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for Negative Binomial(1.3887) family taken to be 1)

    Null deviance: 428.09  on 383  degrees of freedom
Residual deviance: 422.39  on 382  degrees of freedom
AIC: 2342.7

Number of Fisher Scoring iterations: 1


              Theta:  1.389 
          Std. Err.:  0.118 

 2 x log-likelihood:  -2336.683 
> 
> NB1<-glm.nb(LSUM ~ Sex + Treatment)
> summary(NB1)

Call:
glm.nb(formula = LSUM ~ Sex + Treatment, init.theta = 1.5807999, 
    link = log)

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-2.5250  -0.9313  -0.3582   0.3087   2.7566  

Coefficients:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)  2.93889    0.14419  20.382  < 2e-16 ***
Sex         -0.60724    0.09093  -6.678 2.42e-11 ***
TreatmentX  -0.20168    0.09014  -2.237   0.0253 *  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for Negative Binomial(1.5808) family taken to be 1)

    Null deviance: 470.77  on 383  degrees of freedom
Residual deviance: 420.58  on 381  degrees of freedom
AIC: 2303.2

Number of Fisher Scoring iterations: 1


              Theta:  1.581 
          Std. Err.:  0.140 

 2 x log-likelihood:  -2295.168 
> 
> 
> NB2<-glm.nb(LSUM ~ Sex + Treatment + Sex*Treatment)
> summary(NB2)

Call:
glm.nb(formula = LSUM ~ Sex + Treatment + Sex * Treatment, init.theta = 1.583566101, 
    link = log)

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-2.5413  -0.9535  -0.3228   0.2639   2.6825  

Coefficients:
               Estimate Std. Error z value Pr(>|z|)    
(Intercept)      3.0329     0.1915  15.837  < 2e-16 ***
Sex             -0.6735     0.1270  -5.303 1.14e-07 ***
TreatmentX      -0.3959     0.2757  -1.436    0.151    
Sex:TreatmentX   0.1355     0.1818   0.745    0.456    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for Negative Binomial(1.5836) family taken to be 1)

    Null deviance: 471.36  on 383  degrees of freedom
Residual deviance: 420.55  on 380  degrees of freedom
AIC: 2304.6

Number of Fisher Scoring iterations: 1


              Theta:  1.584 
          Std. Err.:  0.140 

 2 x log-likelihood:  -2294.613 
> 