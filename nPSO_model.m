function [x, coords, comm, d] = nPSO_model(N, m, T, gamma, distr, plot_flag)

% Implementation of the Popularity-Similarity-Optimization generative model
% with nonuniform (nPSO) or uniform (PSO) distribution of angular coordinates.

% References:
% 1) nPSO model
%   Muscoloni, A. and Cannistraci, C.V. (2018)
%   "A nonuniform popularity-similarity optimization (nPSO) model
%   to efficiently generate realistic complex networks with communities".
%   New Journal of Physics. doi.org/10.1088/1367-2630/aac06f
% 2) PSO model
%   Papadopoulos, F. et al. (2012)
%   "Popularity versus similarity in growing networks".
%   Nature, 489, 537-540. doi.org/10.1038/nature11459

% Released under MIT License
% Copyright (c) 2017 A. Muscoloni, C. V. Cannistraci

%%% INPUT %%%
% N - number of nodes
% m - half of average node degree (approximately),
%     the number of edges is E = m*(m+1)/2 + (N-m-1)*m
% T - temperature (inversely related to clustering) [T>=0]
% gamma - exponent of the power-law node degree distribution [gamma>=2]
% distr - input defining the distribution of the angular coordinates,
%       it can be:
%       > the value 0 to set a uniform angular distribution
%         (PSO model)
%       > a positive integer indicating the number of components
%         of a Gaussian mixture distribution, in this case the means
%         of the components are equidistant in the angular coordinates space,
%         they have the same standard deviation equal to 1/6 of the distance
%         between two adjacent means, and the mixing proportions
%         are equal for the components
%         (nPSO model with default Gaussian mixture distribution)
%       > a 'gmdistribution' object
%         (nPSO model with custom Gaussian mixture distribution)
%       > a cell having three elements:
%         - vector with evenly spaced points between 0 and 2pi
%         - vector representing the related probability density function
%         - vector with the points representing the center of the communities
%         (nPSO model with custom mixture distribution)
% plot_flag - 1 or 0 to indicate whether the network and the mixture distribution
%             have to be plotted or not (optional, default is 0)
%
%%% OUTPUT %%%
% x - adjacency matrix of the network generated
% coords - polar coordinates (theta, radius) of the network nodes
% comm - community memberships of the network nodes defined as
%        the component of the mixture distribution whose mean
%        is at the lowest angular distance
%        (the output is empty in the PSO model case)
% d - matrix of pairwise hyperbolic distances between the nodes

% NB: for the establishment of the new links in the network we adopt
% the fast implementation referred in the paper as implementation #3.

% check input
narginchk(5,6)
validateattributes(N, {'numeric'}, {'scalar','integer','>=',1,'finite'});
validateattributes(m, {'numeric'}, {'scalar','integer','>=',1,'<',N});
validateattributes(T, {'numeric'}, {'scalar','nonnegative','finite'});
validateattributes(gamma, {'numeric'}, {'scalar','>=',2,'finite'});
if isnumeric(distr)
    validateattributes(distr, {'numeric'}, {'scalar','integer','nonnegative','finite'});
elseif iscell(distr)
    validateattributes(distr, {'cell'}, {'vector','numel',3});
    validateattributes(distr{1}, {'numeric'}, {'vector','>=',0,'<=',2*pi});
    validateattributes(distr{2}, {'numeric'}, {'vector','nonnegative','finite','numel',length(distr{1})});
    validateattributes(distr{3}, {'numeric'}, {'vector','>=',0,'<=',2*pi});
elseif ~isa(distr,'gmdistribution')
    error('Input argument "distr" must be either a numeric scalar, a cell or a gmdistribution');
end
if nargin==5
    plot_flag = 0;
else
    validateattributes(plot_flag, {'numeric'}, {'scalar','integer','>=',0,'<=',1});
end

% initialization
coords = zeros(N,2);
beta = 1 / (gamma - 1);
x = zeros(m*(m+1)/2+(N-m-1)*m,2);
i = 0;

% randomly sample the angular coordinates
% and set the community memberships
if isnumeric(distr) && distr == 0
    % uniform distribution
    uniform = 1;
    coords(:,1) = rand(N,1)*2*pi;
    comm = [];
elseif isnumeric(distr) && distr > 0
    % Gaussian mixture distribution with default settings
    uniform = 0;
    C = distr;
    mu = (0:C-1) .* (2*pi/C);
    sigma = (2*pi/C/6)^2;
    p = repmat(1/C,1,C);
    gmd = gmdistribution(mu', sigma, p);
    coords(:,1) = mod(random(gmd,N),2*pi);
    [~,comm] = min(pi - abs(pi-abs(repmat(coords(:,1),1,C)-repmat(mu,N,1))),[],2);
elseif isa(distr,'gmdistribution')
    % Gaussian mixture distribution with custom settings
    uniform = 0;
    gmd = distr;
    C = gmd.NumComponents;
    mu = gmd.mu';
    coords(:,1) = mod(random(gmd,N),2*pi);
    [~,comm] = min(pi - abs(pi-abs(repmat(coords(:,1),1,C)-repmat(mu,N,1))),[],2);        
elseif iscell(distr)
    % custom mixture distribution
    uniform = 0;
    angles = distr{1};
    angles_prob = distr{2};
    mu = distr{3};
    C = length(mu);
    coords(:,1) = mod(angles(randsample(length(angles),N,true,angles_prob)),2*pi);
    if iscolumn(mu)
        mu = mu';
    end
    [~,comm] = min(pi - abs(pi-abs(repmat(coords(:,1),1,C)-repmat(mu,N,1))),[],2);
end

for t = 2:N
    
    % update the radial coordinates of the existing nodes (popularity fading)
    coords(1:t-1,2) = beta .* (2*log(1:t-1)) + (1-beta)*2*log(t);
    
    % radial coordinate of the new node
    coords(t,2) = 2*log(t);
    
    if t-1 <= m
        % if the existing nodes are equal or lower than the nodes to connect,
        % connect to all of them
        x(i+1:i+t-1,1) = t;
        x(i+1:i+t-1,2) = 1:t-1;
        i = i + t-1;
    else
        % hyperbolic distance of the new node to all the existing nodes
        d = pdist2(coords(t,:), coords(1:t-1,:), @hyperbolic_dist);
        if T == 0
            % connect to the m hyperbolically closest nodes
            [~, idx] = sort(d);
            x(i+1:i+m,1) = t;
            x(i+1:i+m,2) = idx(1:m);
            i = i + m;
        else
            % probability that the new node connects to the existing nodes
            if beta == 1
                Rt = 2*log(t) - 2*log((2*T*log(t))/(sin(T*pi)*m));
            else
                Rt = 2*log(t) - 2*log((2*T*(1 - exp(-(1 - beta)*log(t))))/(sin(T*pi)*m*(1 - beta)));
            end
            p = 1 ./ (1 + exp((d - Rt)./(2*T)));
            
            % nonuniform sampling of m targets according to the probabilities
            idx = datasample(1:t-1, m, 'Replace', false, 'Weights', p);
            x(i+1:i+m,1) = t;
            x(i+1:i+m,2) = idx;
            i = i + m;
        end
    end
end

% create the adjacency matrix from the edge list
x = sparse([x(:,1);x(:,2)],[x(:,2);x(:,1)],1,N,N);

% matrix of pairwise hyperbolic distances between the nodes
if nargout == 4
    d = squareform(pdist(coords, @hyperbolic_dist));
end

if plot_flag
    % plot mixture distribution
    if ~uniform
        figure('color','white')
        if ~iscell(distr)
            [angles, angles_prob] = compute_points_gmd(gmd);
        end
        plot(angles, angles_prob, 'b', 'LineWidth', 2);
        xlim([0,2*pi])
        title('Mixture distribution of angular coordinates', 'FontSize', 14)
    end
    
    % plot network
    figure('color','white')
    plot_hyperbolic_network(x, coords, comm);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_hyperbolic_network(x, coords, comm)

if isempty(comm)
    colours = compute_similarity_colours(coords(:,1));
else
    colours = hsv(max(comm));
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

xlim([-radius, radius])
ylim([-radius, radius])

axis square
axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = hyperbolic_dist(XI, XJ)

% Compute the hyperbolic distance between point XI (1 x 2) and the m points
% stored in XJ (m x 2). The coordinates of these points are in polar format
% (angular coord, radial coord). The resulting similarities are stored in d.
%
% INPUT
% XI -> The polar coordinates of a single point in the Poincaré disk.
% XJ -> The polar coordinates of m points in the Poincaré disk.
%
% OUTPUT
% d -> The hyperbolic distance between point XI and the other m points
%      stored in XJ. The hyperbolic distance between points (Ti, Ri) and
%      (Tj, Rj) in the hyperbolic space H^2 of curvature K = -1,
%      represented by the Poincaré disk is:
%      Dij = arccosh(cosh(Ri)*cosh(Rj) - sinh(Ri)*sinh(Rj)*cos(Tij));
%      where the angular distance is:
%      Tij = pi - |pi - |Ti - Tj||
%
% Copyright (C) Gregorio Alanis-Lobato, 2014
% adapted by Alessandro Muscoloni, 2018

A =  pi - abs(pi - abs(XI(1) - XJ(:,1)));
d = acosh(cosh(XI(2)).*cosh(XJ(:,2)) - sinh(XI(2)).*sinh(XJ(:,2)).*cos(A));
d(isinf(d)) = 0;

% due to numerical approximations, points with a tiny or zero angular
% distance and close radial coordinates could produce a wrong complex
% hyperbolic distance with zero real part. These distances are replaced
% by the radial separation as expected by the theoretical formula.
if ~isreal(d)
    d(imag(d)~=0) = abs(XI(2)-XJ(imag(d)~=0,2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, y] = compute_points_gmd(gmd)

n = 1000;
e = 0.0001;
x = linspace(0,2*pi,n+1);
x = x(1:end-1);
y = pdf(gmd,x');

left = y(1);
right = y(end);
j = 1;
while left > e
    x2 = linspace(-j*2*pi,-(j-1)*2*pi,n+1);
    x2 = x2(1:end-1);
    y2 = pdf(gmd,x2');
    y = y + y2;
    left = y2(1);
    j = j + 1;
end
j = 1;
while right > e
    x2 = linspace(j*2*pi,(j+1)*2*pi,n+1);
    x2 = x2(1:end-1);
    y2 = pdf(gmd,x2');
    y = y + y2;
    right = y2(end);
    j = j + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function colours = compute_similarity_colours(angles)

P = round(2*pi*1000);
angles = mod(round(angles*1000),P)+1;
colours = hsv(P);
colours = colours(angles,:);
