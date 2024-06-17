# Bell-McAffrey standard errors run and are greater than linearized

    Code
      summary(svyglm(api00 ~ api99 + stype, design = dclus1, std.errors = "Bell-McAffrey"))
    Output
      
      Call:
      svyglm(formula = api00 ~ api99 + stype, design = dclus1, std.errors = "Bell-McAffrey")
      
      Survey design:
      dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
      
      Coefficients:
                   Estimate Std. Error t value Pr(>|t|)    
      (Intercept)  99.85645   18.02030   5.541 0.000175 ***
      api99         0.90329    0.02734  33.039 2.33e-12 ***
      stypeH      -19.38726    5.43114  -3.570 0.004398 ** 
      stypeM      -18.15821    6.07011  -2.991 0.012267 *  
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      (Dispersion parameter for gaussian family taken to be 710.3237)
      
      Number of Fisher Scoring iterations: 2
      

# Bell-McAffrey degrees of freedom go down

    Code
      summary(svyglm(api00 ~ api99 + stype, design = dclus1, std.errors = "Bell-McAffrey",
      degf = TRUE))
    Output
      
      Call:
      svyglm(formula = api00 ~ api99 + stype, design = dclus1, std.errors = "Bell-McAffrey", 
          degf = TRUE)
      
      Survey design:
      dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
      
      Coefficients:
                   Estimate Std. Error t value Pr(>|t|)    
      (Intercept)  99.85645   18.02030   5.541  0.00218 ** 
      api99         0.90329    0.02734  33.039 2.39e-07 ***
      stypeH      -19.38726    5.43114  -3.570  0.01452 *  
      stypeM      -18.15821    6.07011  -2.991  0.02824 *  
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      (Dispersion parameter for gaussian family taken to be 710.3237)
      
      Number of Fisher Scoring iterations: 2
      

---

    Code
      summary(svyglm(api00 ~ api99 + stype, design = dclus1, std.errors = "Bell-McAffrey-2",
      degf = TRUE))
    Output
      
      Call:
      svyglm(formula = api00 ~ api99 + stype, design = dclus1, std.errors = "Bell-McAffrey-2", 
          degf = TRUE)
      
      Survey design:
      dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
      
      Coefficients:
                   Estimate Std. Error t value Pr(>|t|)    
      (Intercept)  99.85645   18.02030   5.541  0.00989 ** 
      api99         0.90329    0.02734  33.039  3.8e-05 ***
      stypeH      -19.38726    5.43114  -3.570  0.03412 *  
      stypeM      -18.15821    6.07011  -2.991  0.05386 .  
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      (Dispersion parameter for gaussian family taken to be 710.3237)
      
      Number of Fisher Scoring iterations: 2
      

