function mixture_pdf = create_mixture_gaussian_gamma_pdf(C)

% Generation of Gaussian and Gamma mixture distribution with C components.

% The procedure used is described in Supplementary Information of the paper:
% Muscoloni, A. and Cannistraci, C.V. (2018)
% "A nonuniform popularity-similarity optimization (nPSO) model
% to efficiently generate realistic complex networks with communities".
% New Journal of Physics. doi.org/10.1088/1367-2630/aac06f

%%% INPUT %%%
% C - number of components
%
%%% OUTPUT %%%
% mixture_pdf - a cell having three elements:
%   (1) vector with evenly spaced points between 0 and 2pi
%   (2) vector representing the related probability density function
%   (3) vector with the points representing the center of the communities

% Released under MIT License
% Copyright (c) 2017 A. Muscoloni, C. V. Cannistraci

% set parameters
mu = (0:C-1) .* (2*pi/C);
sigma = 2*pi/C/6;
a = 2;
b = 1/C;

% random ordering of Gaussian and Gamma components (1 = gaussian; 0 = gamma)
mixture = [ones(ceil(C/2),1); zeros(ceil(C/2),1)];
mixture = mixture(randperm(length(mixture)));
mixture = mixture(1:C);

n = 10000;  % number of evenly spaced points between 0 and 2pi
e = 0.0001; % threshold used to truncate the pdf
x = linspace(0,2*pi,n+1);
x = x(1:end-1);
y = zeros(size(x));

for i = 1:C
    
    % make distribution
    if mixture(i)==1
        p = makedist('Normal', mu(i), sigma);
    else
        p = makedist('Gamma', a, b);
    end
    yi = pdf(p,x);
    
    % sum up also the values out of the range [0,2pi]
    % until they are higher than the threshold 'e'
    left = yi(1);
    right = yi(end);
    j = 1;
    while left > e
        x2 = linspace(-j*2*pi,-(j-1)*2*pi,n+1);
        x2 = x2(1:end-1);
        y2 = pdf(p,x2);
        yi = yi + y2;
        left = y2(1);
        j = j + 1;
    end
    j = 1;
    while right > e
        x2 = linspace(j*2*pi,(j+1)*2*pi,n+1);
        x2 = x2(1:end-1);
        y2 = pdf(p,x2);
        yi = yi + y2;
        right = y2(end);
        j = j + 1;
    end
    
    % adjustments of gamma distribution
    if mixture(i)==0
        % translation toward the desired mean
        [~,xmu1] = min(abs(a*b-x));
        [~,xmu2] = min(abs(mu(i)-x));
        d = xmu2 - xmu1;
        if d > 0
            temp = yi;
            yi(1:d) = temp(n-d+1:n);
            yi(d+1:xmu2) = temp(1:xmu1);
            yi(xmu2+1:n) = temp(xmu1+1:n-d);
        elseif d < 0
            d = abs(d);
            temp = yi;
            yi(1:xmu2) = temp(d+1:xmu1);
            yi(xmu2+1:n-d) = temp(xmu1+1:n);
            yi(n-d+1:n) = temp(1:d);
        end
        % reflection with respect to the mean (random option)
        if rand > 0.5
           if xmu2 > n/2
               temp = yi;
               yi(xmu2:-1:2*xmu2-n) = temp(xmu2:n);
               yi(xmu2:n) = temp(xmu2:-1:2*xmu2-n);
               yi(1:2*xmu2-n-1) = temp(2*xmu2-n-1:-1:1);
           else
               temp = yi;
               yi(2*xmu2-1:-1:xmu2) = temp(1:xmu2);
               yi(xmu2:-1:1) = temp(xmu2:2*xmu2-1);
               yi(2*xmu2:n) = temp(n:-1:2*xmu2);
           end
        end
    end 
    y = y + yi;
end
y = y ./ C;

mixture_pdf = {x; y; mu};
