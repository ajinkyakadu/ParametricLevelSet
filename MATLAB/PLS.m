function [f, g, m] = PLS(x, funObj, A, options)
% Generates the Parametric Level-Set function for any function of R^n. 
% The optimization problem is reformulated as finding x in R^m so that the 
% objective defined by the function F is minimized. The model vector on R^n 
% is represented by a vector of RBFs included in the kernel matrix A.
%
% INPUTS:
%   x - vector defined on R^m, representing the Parametric Level-Set.
%   funObj - function handle for solving the inverse problem.
%   A - RBF kernel matrix with dimensions [n x m].
%   options - struct containing parameters with the following possible fields:
%      m0 - background model vector of size n (default: zeros(n,1)).
%      m1 - parameter for the level-set function (default: 1).
%      kappa - parameter controlling the width of the Heaviside function 
%           (default: 0.1).
%      hopt - struct containing parameters for the Heaviside function.
%
% OUTPUTS:
%   f - function value from the Parametric Level-Set at x.
%   g - gradient vector of size m at x.
%   m - resulting model vector.
%
% Supported Input options:
%   hopt - struct containing parameters for the Heaviside function with the
%       following possible fields:
%       epsi - parameter determining the width of the Heaviside function 
%           (default: 0.1).
%
% Author: Ajinkya Kadu
%         Mathematical Institute,
%         Utrecht University, The Netherlands
%
% Date : September 2016

if nargin < 4
    options = [];
end

% Set default values for non-existent or blank fields in the options struct
m0    = GetOptions(options, 'm0', zeros(size(A,1),1));
m1    = GetOptions(options, 'm1', 1);
kappa = GetOptions(options, 'kappa', 0.1);

Ax = A*x;

% Calculate the width of the Heaviside function
options.hopt.epsi = kappa * 0.5 * (max(Ax) - min(Ax));
[h,d] = HeavisideFunction(Ax,options.hopt);

% Calculate the resulting model vector
m = m0 .* (1 - h) + h * m1;

% Calculate the function value and gradient vector
[f,g0]  = funObj(m);
g = A' * (g0 .* d .* (m1 - m0));

end