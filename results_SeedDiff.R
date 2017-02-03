kruskal.test(SeedDiff~Treatment, data=gPaired2)

	Kruskal-Wallis rank sum test

data:  SeedDiff by Treatment
Kruskal-Wallis chi-squared = 1.3179, df = 1, p-value = 0.251


> GLM1<-glm(SeedDiff~Treatment, data=gPaired2)
> summary(GLM1)

Call:
glm(formula = SeedDiff ~ Treatment, data = gPaired2)

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-3.4566  -0.4010  -0.1931   0.4568   2.5028  

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)   0.1939     0.1319   1.470    0.144
TreatmentX    0.1543     0.1797   0.859    0.392

(Dispersion parameter for gaussian family taken to be 1.026597)

    Null deviance: 130.11  on 127  degrees of freedom
Residual deviance: 129.35  on 126  degrees of freedom
AIC: 370.59

Number of Fisher Scoring iterations: 2
