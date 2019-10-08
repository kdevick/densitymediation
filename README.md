************************************
# Density Mediation 



Code developed by:   Katrina Devick (devick.katrina@mayo.edu)  
Last updated:   7 October 2019



******************
## Getting Started


The files density_mediation.R and density_mediation_source.R that can be downloaded from the repository implement our proposed density mediation approach from ArXiv article ["The Role of Body Mass Index at Diagnosis on Black-White Disparities in Colorectal Cancer Survival: A Density Regression Mediation Approach"](https://arxiv.org/abs/1812.02829). The file code_to_extract_cond_densities.R demonstrates how to obtain conditional density and CDFs for each MCMC iteration from a LDDPdensity fit. 



*******************
## Important Note


Prior to implementing the code, you must download the file LDDPdensity.f, replace the current LDDPdensity.f file in the DPpackage with this file, and then recompile and load this new DPpackage. You can easily do this using devtools and downloading the package "DPpackage-modified" in this repository with the following code: 


```{r, eval=FALSE}
install.packages("devtools")
library(devtools)

install_github("kdevick/densitymediation/DPpackage-modified")
library(DPpackage)
```

Without this step, you are unable to extract the conditional density and CDFs for each MCMC iteration, which are needed to obtain posterior samples of the residual disparity (RD). Please note, if you are using a Mac, you may get the following error when executing the above code:

```{r, eval=FALSE}
make: gfortran: No such file or directory
make: *** [BDPdensity.o] Error 1
ERROR: compilation failed for package ‘DPpackage’
```

This error is because newer versions of R for a Mac no longer include the gfortran tools needed to compile a R package from the source code, if the source code is written in Fortran. In order to compile this modified DPpackage with a newer version of R on a Mac, you need to first download the Fortran binary files found [here](https://github.com/fxcoudert/gfortran-for-macOS/releases), and then restart your device. 



