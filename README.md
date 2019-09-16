# L1MLR
## L1-penalized multinomial logistic regression models for unmatched case control studies

Methods MultinomLogist_StdLasso and MultinomLogist_SymLasso, developped by Ballout, Garcia and Viallon *, are used to estimate L1-penalized multinomial logistic regression models for unmatched case control studies. Functions corresponding to these methods are available in MultinomLogist_StdLasso.R and MultinomLogist_SymLasso respectively).

**MultinomLogist_StdLasso** returns estimates using the standard parametrization based on a a priori selection of a reference category.

**MultinomLogist_SymLasso** returns estimates using the symmetric formulation, as implemented in the glmnet package. This is equivalent to using the standard formulation and applying a data shared lasso penalty.

See Illustr.R for a simple example illustrating the use of these functions.

See (https://arxiv.org/abs/1901.01583) [1] for more details.
## Packages required 


```
glmnet, Matrix and Rcpp.
```



## Usage
##### MultinomLogist_RefLasso and MultinomLogist_SymLasso
#### Arguments
* **X**        : Input matrix, of dimension n x p, where n is the number of observations and p is the number of variables; each row is an observation vector.  
* **Y**        : Numerical or factorial vector at K different values. For the MultinomLogist_RefLasso method, the reference is chosen as the first level (i.e. levels (factor (Y)) [1]).
* **intercept**    : Logical. Should an intercept be included in the model - default is FALSE.  
* **method**        : Character string, specifies the tuning parameter selection method to be used. Choices are "BIC", "BIC-H" and/or "CV".  
"BIC" :  specifies the **B**ayesian **I**nformation **C**riterion;  
"BIC-H":  specifies the **B**ayesian **I**nformation **C**riterion, adapting the Hybrid-OLS idea (Efron and others 2004) [2];  
"CV"  :  specifies the **C**ross **V**alidation technique;  

* **version_cpp**      : This is for the "CondLogist_RefLasso" method. If TRUE, performs the iterative algorithm with a function using the Cpp language. This makes the function faster - default is FALSE.

#### Value
* **BIC**         : Matrix of parameters obtained for the lambda value that minimizes the BIC.    
* **BIC_H**       : Matrix of parameters obtained for the lambda value that minimizes the BIC-H.  .   
* **CV**       : Matrix of parameters obtained for the lambda value that minimizes the mean cross validation error.   
   
**ALL outputs are matrices, of dimension p x (K-1), where p is the number of variables and K is the number of levels of Y. The K-1 columns of each matrix corresponds to vector of parameters (log odds ratios) comparing each category to the K-th one.**  




## References

[1] Ballout, Nadim, Garcia, Cedric, et Viallon, Vivian. Sparse estimation for case-control studies with multiple subtypes of cases. arXiv preprint arXiv:1901.01583, 2019.

[2] Efron, Bradley; Hastie, Trevor; Johnstone, Iain; Tibshirani, Robert. Least angle regression. Ann. Statist. 32 (2004), no. 2, 407--499.

