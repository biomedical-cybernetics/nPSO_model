function [x, coords, comm] = nPSO_model(N, m, T, gamma, gmd, method, plot_options)

% Implementation of the Popularity-Similarity-Optimization generative model
% with uniform (PSO) or non-uniform (nPSO) distribution of angular coordinates.

% References:
% 1) PSO model
%   Papadopoulos, F. et al. (2012)
%   "Popularity versus similarity in growing networks". Nature, 489, 537-540.
% 2) nPSO model
%   Muscoloni, A. and Cannistraci, C.V. (2017)
%   "A nonuniform popularity-similarity optimization (nPSO) model
%   to efficiently generate realistic complex networks with communities".
%   arXiv:1707.07325

% Authors:
% The implementation 1 of the PSO model has been coded by Gregorio Alanis-Lobato (2014).
% The variant with non-uniform distribution of angular coordinates,
% as well as the implementations 2 and 3, have been coded by Alessandro Muscoloni (2017).

% Released under MIT License
% Copyright (c) 2017 A. Muscoloni, C. V. Cannistraci

%%% INPUT %%%
% N - number of nodes
% m - half of average degree (defines the number of edges E = N*m)
% T - temperature (inversely related to clustering)
% gamma - exponent of the power-law node degree distribution
% gmd - input defining the distribution of the angular coordinates;
%       it can be:
%       1) the value 0 to set a uniform angular distribution
%       2) a positive integer indicating the number of components
%          of a Gaussian mixture distribution, in this case the means
%          of the components are equidistant in the angular coordinates space,
%          they have the same standard deviation equal to 1/6
%          of the distance between two adjacent means,
%          and the mixing proportions are equal for the components.
%       3) a matrix having as many rows as the number of components
%          of a Gaussian mixture distribution and three columns
%          indicating for each component respectively
%          > the mean of the component (between 0 and 2pi)
%          > the variance of the component
%            if at least one variance is set to 0 the default values are used,
%            with all the components having the same standard deviation,
%            equal to 1/6 of the distance between two adjacent means
%            considering as if the means were equidistant
%          > the mixing proportion of the component
%            if at least one proportion is set to 0 the default values are used,
%            with all the components having the same mixing proportion
%          for more details refer to:
%          https://mathworks.com/help/stats/gmdistribution.html
% method - number indicating which implementation to use for creating
%          the connections when a new node enters in the network:
%          1 - original implementation
%          2 - original implementation with speedup v1
%          3 - original implementation with speedup v2 (recommended)
% plot_options - vector of three columns having as values 1 or 0 to indicate respectively:
%          1) if the Gaussian mixture distribution has to be plotted or not
%          2) if the network has to be plotted or not
%          3) if also an equivalent network with uniform angular distribution
%             has to be plotted or not for comparison
%
% Note: the function can be called
% - with 7 input arguments
%   nPSO_model(N, m, T, gamma, gmd, method, plot_options)
% - with 6 input arguments
%   nPSO_model(N, m, T, gamma, gmd, method)
%   in this case it plots only the network
% - with 5 input arguments
%   nPSO_model(N, m, T, gamma, gmd)
%   in this case it plots only the network
%   and the implementation 3 is used
%
%%% OUTPUT %%%
% x - adjacency matrix of the generated network
% coords - polar coordinates (theta, radius) of the network nodes
% comm - community memberships of the network nodes
%        defined as the closest component of the Gaussian mixture
%        distribution (the output is all-zero in the uniform case)

% check input
if nargin < 5 || nargin > 7
   error(['Input arguments must be between 5 and 7.' char(10) ...
      'Possible usages:' char(10) ...
      'nPSO_model(N, m, T, gamma, gmd, method, plot_options)' char(10) ...
      'nPSO_model(N, m, T, gamma, gmd, method)'  char(10) ...
      'nPSO_model(N, m, T, gamma, gmd)']);
end
validateattributes(N, {'numeric'}, {'scalar','integer','>=',1,'finite'});
validateattributes(m, {'numeric'}, {'scalar','integer','>=',1,'<=',N});
validateattributes(T, {'numeric'}, {'scalar','nonnegative','finite'});
validateattributes(gamma, {'numeric'}, {'scalar','>',2,'finite'});
if isscalar(gmd)
    validateattributes(gmd, {'numeric'}, {'scalar','integer','nonnegative','finite'});
elseif isnumeric(gmd)
    validateattributes(gmd, {'numeric'}, {'ndims',2,'ncols',3});
    validateattributes(gmd(:,1), {'numeric'}, {'>=',0,'<=',2*pi});
    validateattributes(gmd(:,2), {'numeric'}, {'nonnegative','finite'});
    validateattributes(gmd(:,3), {'numeric'}, {'nonnegative','finite'});
else
    error('Input argument "gmd" must be either a scalar or a 2D matrix');
end
if ~exist('method','var')
    method = 3;
else
    validateattributes(method, {'numeric'}, {'scalar','integer','>=',1,'<=',3});
end    
if ~exist('plot_options', 'var')
    plot_gmd = 0;
    plot_network = 1;
    plot_uniform = 0;
else
    validateattributes(plot_options, {'numeric'}, {'size',[1,3],'integer','>=',0,'<=',1});
    plot_gmd = plot_options(1);
    plot_network = plot_options(2);
    plot_uniform = plot_options(3);
end

% initialization
x = sparse(N,N);
coords = zeros(N,2);
comm = zeros(N,1);
beta = 1 / (gamma - 1);

if isscalar(gmd) && gmd == 0
    % uniform angular distribution
    uniform = 1;
elseif isscalar(gmd) && gmd > 0
    % default Gaussian mixture angular distribution
    k = gmd;
    uniform = 0;
    mu = (0:k-1) .* (2*pi/k);
    sigma = (2*pi/k/6)^2;
    p = repmat(1/k,1,k);
    gmd = gmdistribution(mu', sigma, p);
else
    % custom Gaussian mixture angular distribution
    uniform = 0;
    k = size(gmd,1);
    if any(gmd(:,2)==0)
        gmd(:,2) = (2*pi/k/6)^2;
    end
    if any(gmd(:,3)==0)
        gmd(:,3) = repmat(1/k,1,k);
    end
    gmd = gmdistribution(gmd(:,1),reshape(gmd(:,2),1,1,k),gmd(:,3)');
end

for t = 1:N
    
    % move the existing nodes to their new radial coordinates
    % to simulate popularity fading
    coords(1:(t-1),2) = beta .* (2*log(1:(t-1))) + (1-beta)*2*log(t);
    
    % introduce the new node
    % radial coordinate
    coords(t,2) = 2*log(t);
    % angular coordinate
    if uniform
        coords(t,1) = rand(1)*2*pi;
    else
        y = random(gmd);
        if y < 0
            coords(t,1) = 2*pi - mod(abs(y),2*pi);
        elseif y > 2*pi
            coords(t,1) = mod(y,2*pi);
        else
            coords(t,1) = y;
        end
        
        % assign community membership
        [~,comm(t)] = min(pi - abs(pi-abs(coords(t,1)-gmd.mu)));
     end
    
    % hyperbolic distance of the new node to all existing nodes
    d = pdist2(coords(t,:), coords(1:(t-1),:), @hyperbolic_dist);
    
    if T > 0 && t > 1
        if (t-1) <= m
            % if the existing nodes are equal or lower than the nodes to connect,
            % connect all of them
            x(t,1:(t-1)) = 1;
        else
            % probability that the new node connects to the existing nodes
            Rt = 2*log(t) - 2*log((2*T*(1 - exp(-(1 - beta)*log(t))))/(sin(T*pi)*m*(1 - beta)));
            p = 1 ./ (1 + exp((d - Rt)./(2*T)));
            targets = 1:(t-1);
            m_count = 0;
            while m_count < m
                if method == 1
                    rnd_p = rand(1);
                    idx = randi(length(targets));
                    rnd_node = targets(idx);
                    if p(rnd_node) > rnd_p
                        x(t,rnd_node) = 1;
                        m_count = m_count + 1;
                        targets(idx) = [];
                    end
                elseif method == 2
                    if length(targets) > m
                        idx = randsample(length(targets),m);
                        rnd_p = rand(1,m) * max(p(targets));
                    else
                        idx = 1:length(targets);
                        rnd_p = rand(1,length(targets)) * max(p(targets));
                    end
                    idx = idx(p(targets(idx)) > rnd_p);
                    if length(idx) > m - m_count
                        idx = randsample(idx,m - m_count);
                    end
                    x(t,targets(idx)) = 1;
                    targets(idx) = [];
                    m_count = m_count + length(idx);
                elseif method == 3
                    idx = unique(randsample(length(targets),m,1,p(targets)));
                    if length(idx) > m - m_count
                        idx = randsample(idx,m - m_count);
                    end
                    x(t,targets(idx)) = 1;
                    targets(idx) = [];
                    m_count = m_count + length(idx);
                end
            end
        end
    else
        % if T is 0 connect to the m hyperbolically closest nodes
        [~, idx] = sort(d);
        if numel(idx) < m
            x(t,idx) = 1;
        else
            x(t,idx(1:m)) = 1;
        end
    end
end

x = max(x, x');

% plot Gaussian mixture distribution
if ~uniform && plot_gmd
    P = 100;
    X2 = linspace(0,2*pi,P);
    X1 = X2 - 2*pi;
    X3 = X2 + 2*pi;
    X_gmd = [X1(1:end-1) X2(1:end-1) X3(1:end-1)]; clear X1 X2 X3;
    Y_gmd = pdf(gmd, X_gmd');
    Y_gmd = Y_gmd(1:P-1) + Y_gmd(P:2*(P-1)) + Y_gmd(2*(P-1)+1:3*(P-1));
    
    figure('color','white')
    plot(X_gmd(P:2*(P-1)),Y_gmd, 'b', 'LineWidth', 2);
    xlim([X_gmd(P-1),X_gmd(2*(P-1)+1)])
    title('Gaussian mixture distribution of angular coordinates', 'FontSize', 14)
end

% plot network
if plot_network
    figure;
    set(gcf,'color','white')

    % compare network with uniform distribution of angular coordinates
    if ~uniform && plot_uniform
        subplot(1,2,1)
        plot_hyperbolic_network(x, coords, comm, 'non-uniform PSO network');

        subplot(1,2,2)
        [x_uni, coords_uni, comm_uni] = nPSO_model(N, m, T, gamma, 0, method, [0 0 0]);
        plot_hyperbolic_network(x_uni, coords_uni, comm_uni, 'uniform PSO network');
    else
        plot_hyperbolic_network(x, coords, comm, 'PSO network');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_hyperbolic_network(x, coords, comm, string_title)

C = length(unique(comm));
if C == 1
    colormap('hsv')
    colours = coords(:,1);
else
    colours = hsv(C);
    colours = colours(comm,:);
end

% plot circle
hold on
radius = max(coords(:,2));
viscircles([0 0], radius, 'EdgeColor', 'k');

[cart_coords(:,1),cart_coords(:,2)] = pol2cart(coords(:,1),coords(:,2));

% plot links
[h1,h2] = gplot(x, cart_coords,'k');
plot(h1, h2, 'Color',[0.80,0.80,0.80], 'LineWidth',1)

% plot nodes
hold on
scatter(cart_coords(:,1),cart_coords(:,2),200,colours,'filled','MarkerEdgeColor','k');

xlim([-11/10*radius, 11/10*radius])
ylim([-11/10*radius, 11/10*radius])

title(string_title, 'FontSize', 14)

axis square
axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = hyperbolic_dist(XI, XJ)
%Computes the hyperbolic distance between point XI (1 x 2) and the m points
%stored in XJ (m x 2). The coordinates of these points are polar in the
%format (angular coord, radial coord). The resulting similarities are
%stored in d.
%
%INPUT
%   XI -> The polar coordinates of a single point in the Poincaré disc.
%   XJ -> The polar coordinates of m points in the Poincaré disc.
%
%OUTPUT
%   d -> The hyperbolic distance between point XI and the other m points
%        stored in XJ. The hyperbolic distance between points (Ti, Ri) and
%        (Tj, Rj) in the hyperbolic space H^2 of curvature K = -1,
%        represented by the Poincaré disc is:
%
% Dij = arccosh(cosh(Ri)*cosh(Rj) - sinh(Ri)*sinh(Rj)*cos(Tij));
%
%        with Tij = pi - |pi - |Ti - Tj||
%
% Copyright (C) Gregorio Alanis-Lobato, 2014

A =  pi - abs(pi - abs(XI(1) - XJ(:,1))); %angular separation between points
d = acosh(cosh(XI(2)).*cosh(XJ(:,2)) - sinh(XI(2)).*sinh(XJ(:,2)).*cos(A));
d(isinf(d)) = 0;

% due to numerical approximations, points with a tiny or zero angular
% separation and close radial coordinates could produce a wrong complex
% hyperbolic distance with zero real part. These distances are replaced
% by the radial separation as expected by the theoretical formula.
if ~isreal(d)
   d(imag(d)~=0) = abs(XI(2)-XJ(imag(d)~=0,2)); 
end
