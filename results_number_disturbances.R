> #============================================================================
> #      non parameteric test of disturbance by treatment
> #           -Total disturbances ~ Treatment
> #           -Difference between shaded and unshaded disturbances by treatment
> #============================================================================
> 
> kruskal.test(TotalDist~Treatment, data=gd4Paired)

	Kruskal-Wallis rank sum test

data:  TotalDist by Treatment
Kruskal-Wallis chi-squared = 3.5127, df = 1, p-value = 0.0609

> 
> TotalDist_S_minus_u<-gd4Paired$TotalDistS-gd4Paired$TotalDistU
> 
> kruskal.test(TotalDist_S_minus_u~gd4Paired$Treatment)

	Kruskal-Wallis rank sum test

data:  TotalDist_S_minus_u by gd4Paired$Treatment
Kruskal-Wallis chi-squared = 0.18043, df = 1, p-value = 0.671

> 
> #============================================================================
> #      GLMs
> #       -Poisson model of counts of total disturbances by treatment
> #       -same with fixed effect for Site
> #============================================================================
> 
> GLM3<- glm(TotalDistS~Treatment, data=gd4Paired, family = "poisson")
> summary(GLM3)

Call:
glm(formula = TotalDistS ~ Treatment, family = "poisson", data = gd4Paired)

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-2.3152  -0.8098   0.1918   0.6264   1.2098  

Coefficients:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)   0.7138     0.1429   4.996 5.84e-07 ***
TreatmentX    0.2721     0.1880   1.447    0.148    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for poisson family taken to be 1)

    Null deviance: 37.315  on 48  degrees of freedom
Residual deviance: 35.197  on 47  degrees of freedom
AIC: 163.87

Number of Fisher Scoring iterations: 5

> GLM4<- glm(TotalDistS~Treatment + Site, data=gd4Paired, family = "poisson")
> summary(GLM4)

Call:
glm(formula = TotalDistS ~ Treatment + Site, family = "poisson", 
    data = gd4Paired)

Deviance Residuals: 
     Min        1Q    Median        3Q       Max  
-2.46769  -0.44809   0.01409   0.45638   0.99149  

Coefficients:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)  0.84598    0.18602   4.548 5.42e-06 ***
TreatmentX   0.26744    0.18805   1.422   0.1550    
SiteH       -0.02295    0.21332  -0.108   0.9143    
SiteT       -0.42900    0.24284  -1.767   0.0773 .  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for poisson family taken to be 1)

    Null deviance: 37.315  on 48  degrees of freedom
Residual deviance: 31.232  on 45  degrees of freedom
AIC: 163.91

Number of Fisher Scoring iterations: 5

> 