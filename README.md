************************************
# Density Mediation 



Code developed by:   Katrina Devick (kdevick@hsph.harvard.edu)  
Last updated:   10 June 2019



******************
## Getting Started


The files density_mediation.R and density_mediation_source.R that can be downloaded from the repository implement our proposed density mediation approach from ArXiv article ["The Role of Body Mass Index at Diagnosis on Black-White Disparities in Colorectal Cancer Survival: A Density Regression Mediation Approach"](https://arxiv.org/abs/1812.02829). The file code_to_extract_cond_densities.R demonstrates how to obtain conditional density and CDFs for each MCMC iteration from a LDDPdensity fit. 



*******************
## Important Note


Prior to implementing the code, you must download the file LDDPdensity.f, replace the current LDDPdensity.f file in the DPpackage with this file, and then recompile and load this new DPpackage. Without this step, you are unable to extract the conditional density and CDFs for each MCMC iteration, which are needed to obtain posterior samples of the residual disparity (RD). 


