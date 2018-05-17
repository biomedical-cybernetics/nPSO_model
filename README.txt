################################################################################

> nPSO_model.m

Implementation of the Popularity-Similarity-Optimization generative model
with nonuniform (nPSO) or uniform (PSO) distribution of angular coordinates.

### REFERENCES ###

1) nPSO model
Muscoloni, A. and Cannistraci, C.V. (2018)
"A nonuniform popularity-similarity optimization (nPSO) model
to efficiently generate realistic complex networks with communities".
New Journal of Physics. doi.org/10.1088/1367-2630/aac06f

2) PSO model
Papadopoulos, F. et al. (2012)
"Popularity versus similarity in growing networks".
Nature, 489, 537-540. doi.org/10.1038/nature11459

Released under MIT License
Copyright (c) 2017 A. Muscoloni, C. V. Cannistraci

Usage:
[x, coords, comm] = nPSO_model(N, m, T, gamma, distr, plot_flag)

### INPUT ###
N - number of nodes
m - half of average node degree (approximately),
    the number of edges is E = m*(m+1)/2 + (N-m-1)*m
T - temperature (inversely related to clustering) [T>=0]
gamma - exponent of the power-law node degree distribution [gamma>2]
distr - input defining the distribution of the angular coordinates,
      it can be:
      > the value 0 to set a uniform angular distribution
        (PSO model)
      > a positive integer indicating the number of components
        of a Gaussian mixture distribution, in this case the means
        of the components are equidistant in the angular coordinates space,
        they have the same standard deviation equal to 1/6 of the distance
        between two adjacent means, and the mixing proportions
        are equal for the components
        (nPSO model with default settings)
      > a cell having three elements:
        - vector with evenly spaced points between 0 and 2pi
        - vector representing the related probability density function
        - vector with the points representing the center of the communities
        (nPSO model with custom settings)
plot_flag - 1 or 0 to indicate whether the network and the mixture distribution
            have to be plotted or not (optional, default is 0)

### OUTPUT ###
x - adjacency matrix of the network generated
coords - polar coordinates (theta, radius) of the network nodes
comm - community memberships of the network nodes defined as
       the component of the mixture distribution whose mean
       is at the lowest angular distance
       (the output is empty in the PSO model case)

NB: for the establishment of the new links in the network we adopt
the fast implementation referred in the paper as implementation #3.	   

################################################################################

> create_mixture_gaussian_gamma_pdf.m

Generation of Gaussian and Gamma mixture distribution with C components.

### REFERENCES ###

The procedure used is described in Supplementary Information of the paper:
Muscoloni, A. and Cannistraci, C.V. (2018)
"A nonuniform popularity-similarity optimization (nPSO) model
to efficiently generate realistic complex networks with communities".
New Journal of Physics. doi.org/10.1088/1367-2630/aac06f

Released under MIT License
Copyright (c) 2017 A. Muscoloni, C. V. Cannistraci

Usage: mixture_pdf = create_mixture_gaussian_gamma_pdf(C)

### INPUT ###
C - number of components

### OUTPUT ###
mixture_pdf - a cell having three elements:
 (1) vector with evenly spaced points between 0 and 2pi
 (2) vector representing the related probability density function
 (3) vector with the points representing the center of the communities

 
################################################################################

> example_usage.m

Example usage of the functions "nPSO_model" and "create_mixture_gaussian_gamma_pdf".
The script generates and plots three examples:
- nPSO model with 4 communities and default settings (Gaussian mixture distribution)
- nPSO model with 8 communities and custom settings (Gaussian and Gamma mixture distribution)
- PSO model

################################################################################

### CONTACT ###

For any problem, please contact:
Alessandro Muscoloni: alessandro.muscoloni@gmail.com
Carlo Vittorio Cannistraci: kalokagathos.agon@gmail.com