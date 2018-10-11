# Nonuniform PSO (nPSO) model

## Reference

A. Muscoloni, and C.V. Cannistraci, "A nonuniform popularity-similarity optimization (nPSO) model to efficiently generate realistic complex networks with communities", New Journal of Physics, 2018.  
https://doi.org/10.1088/1367-2630/aac06f

## Files description

* *nPSO_model.m*  
  Implementation of the Popularity-Similarity-Optimization generative model with nonuniform (nPSO) or uniform (PSO) distribution of angular coordinates.

* *create_mixture_gaussian_gamma_pdf.m*  
  Generation of Gaussian and Gamma mixture distribution with C components, as described in the Supplementary Information of the reference.

* *example_usage.m*  
  Example usage of the previous functions, the script generates and plots four examples:
  - nPSO model with 4 communities and default Gaussian mixture distribution
  - nPSO model with 4 communities and custom Gaussian mixture distribution
  - nPSO model with 8 communities and custom mixture distribution (Gaussian and Gamma mixture distribution)
  - PSO model

## Contact

For any problem, please contact:
* Alessandro Muscoloni: alessandro.muscoloni@gmail.com
* Carlo Vittorio Cannistraci: kalokagathos.agon@gmail.com