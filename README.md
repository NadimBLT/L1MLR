# L1MLR
## L1-penalized multinomial logistic regression models for unmatched case control studies

Methods MultinomLogist_RefLasso and MultinomLogist_SymLasso, developped by Ballout, Garcia and Viallon *, are used to estimate L1-penalized multinomial logistic regression models for unmatched case control studies. Functions corresponding to these methods are available in MultinomLogist_RefLasso.R and MultinomLogist__Ref_SymLasso respectively).

**MultinomLogist_RefLasso** estimates the models in an arbitrary way by fitting L1-penalized multinomial logistic regression after choosing an arbitrary reference category.

**MultinomLogist_SymLasso** estimates the models in a joint way by fitting L1-penalized multinomial logistic regression using the data shared lasso penalty (Viallon *). This method by-passes the arbitrary choice of the reference category.

See Illustr.R for a simple example illustrating the use of these functions.

See https:...... for more details.
## Packages required 


```
glmnet, Matrix and Rcpp.
```



## Usage
##### MultinomLogist_RefLasso and MultinomLogist_SymLasso
#### Arguments
* **X**        : Input matrix, of dimension n x p, where n is the number of observations and p is the number of variables; each row is an observation vector.  
* **Y**        : Numerical or factorial vector at K different values. For the MultinomLogist_RefLasso method, the reference is chosen as the first level (i.e. levels (factor (Y)) [1]).
* **intercept**    : Vector defining the pairs; each pair composed by a case and his matched control.  
* **method**        : Character string, specifies the tuning parameter selection method to be used. Choices are "BIC", "BIC-H" and/or "CV".  
"BIC" :  specifies the **B**ayesian **I**nformation **C**riterion;  
"BIC-H":  specifies the **B**ayesian **I**nformation **C**riterion **H**ybrid;  
"CV"  :  specifies the **C**ross **V**alidation technique;  

* **version_cpp**      : This is for the "CondLogist_RefLasso" method. If TRUE, performs the iterative algorithm with a function using the Cpp language. This makes the function faster - default is FALSE.

#### Value
* **BIC**         : Estimated models returned by the BIC method, models that have the minmum values of BIC.    
* **BIC_H**       : Estimated models returned by the BIC-R method, models that have the minmum values of BIC in Hybrid version.  .   
* **CV**       : Estimated models returned by the CV method, models that have the minmum values of mean cross validation error.   
   
**ALL outputs are matrices, of dimension p x (K-1), where p is the number of variables and K is the number of levels of Y.**  
