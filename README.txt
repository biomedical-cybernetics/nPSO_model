Implementation of the Popularity-Similarity-Optimization generative model
with uniform (PSO) or non-uniform (nPSO) distribution of angular coordinates.


### REFERENCES ###

1) PSO model
Papadopoulos, F. et al. (2012)
"Popularity versus similarity in growing networks". Nature, 489, 537-540.

2) nPSO model
Muscoloni, A. and Cannistraci, C.V. (2017)
"A nonuniform popularity-similarity optimization (nPSO) model
to efficiently generate realistic complex networks with communities".
arXiv:1707.07325


### INPUT ###
N - number of nodes
m - half of average degree (defines the number of edges E = N*m)
T - temperature (inversely related to clustering)
gamma - exponent of the power-law node degree distribution
gmd - input defining the distribution of the angular coordinates;
      it can be:
      1) the value 0 to set a uniform angular distribution
      2) a positive integer indicating the number of components
         of a Gaussian mixture distribution, in this case the means
         of the components are equidistant in the angular coordinates space,
         they have the same standard deviation equal to 1/6
         of the distance between two adjacent means,
         and the mixing proportions are equal for the components.
      3) a matrix having as many rows as the number of components
         of a Gaussian mixture distribution and three columns
         indicating for each component respectively
         > the mean of the component (between 0 and 2pi)
         > the variance of the component
           if at least one variance is set to 0 the default values are used,
           with all the components having the same standard deviation,
           equal to 1/6 of the distance between two adjacent means
           considering as if the means were equidistant
         > the mixing proportion of the component
           if at least one proportion is set to 0 the default values are used,
           with all the components having the same mixing proportion
         for more details refer to:
         https://mathworks.com/help/stats/gmdistribution.html
method - number indicating which implementation to use for creating
         the connections when a new node enters in the network:
         1 - original implementation
         2 - original implementation with speedup v1
         3 - original implementation with speedup v2 (recommended)
plot_options - vector of three columns having as values 1 or 0 to indicate respectively:
         1) if the Gaussian mixture distribution has to be plotted or not
         2) if the network has to be plotted or not
         3) if also an equivalent network with uniform angular distribution
            has to be plotted or not for comparison

Note: the function can be called
- with 7 input arguments
  PSO_model(N, m, T, gamma, gmd, method, plot_options)
- with 6 input arguments
  PSO_model(N, m, T, gamma, gmd, method)
  in this case it plots only the network
- with 5 input arguments
  PSO_model(N, m, T, gamma, gmd)
  in this case it plots only the network
  and the implementation 3 is used

  
### OUTPUT ###
x - adjacency matrix of the generated network
coords - polar coordinates (theta, radius) of the network nodes
comm - community memberships of the network nodes
       defined as the closest component of the Gaussian mixture
       distribution (the output is all-zero in the uniform case)

	   
### CONTACT ###

For any problem, please contact:
Alessandro Muscoloni: alessandro.muscoloni@gmail.com
Carlo Vittorio Cannistraci: kalokagathos.agon@gmail.com