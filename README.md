Ordinal regression analysis for continuous scales
==================

Ordinal regression analysis is a convenient tool for analyzing ordinal response variables in the presence of covariates. We extend this methodology to the case of continuous self-rating scales such as the Visual Analog Scale (VAS) used in pain assessment, or the Linear Analog Self-Assessment (LASA) scales in quality of life studies. Subjects are typically given a linear scale of 100 mm and asked to put a mark where they perceive themselves. These scales  measure subjects' perception of an intangible quantity, and cannot be handled as ratio variables because of their inherent nonlinearity.  Instead we treat them as ordinal variables, measured on a continuous scale. We express  the likelihood in terms of a function (the ``g function'') connecting the scale with an underlying continuous latent  variable. In the current version the g function is taken as the generalized logistic function (Richards 1959). This has 3 parameters: *M*, the offset, *B*, the slope, and *T*, the symmetry of the curve. The link function is the inverse of the CDF of the assumed underlying distribution of the latent variable. Currently the logit link, which corresponds to a standard logistic distribution, is implemented (this implies a proportional odds model). The likelihood is maximized using *optim {stats}* with a quasi-Newton method (*"BFGS"*). Fixed-effects models are implemented in the function *ocm*, and mixed models in  *ocmm*. 

### References

* Manuguerra M, Heller GZ (2010). Ordinal Regression Models for Continuous Scales, *The International Journal of Biostatistics*: 6(1), Article 14.
* Richards, F. (1959). A flexible growth function for empirical use, *Journal of Experimental Botany*, 10, 290-301.
