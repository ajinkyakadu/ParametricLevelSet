function [A,nr] = GenerateKernel(x, z, options)
%GenerateKernel Generates Radial Basis Function (RBF) Kernel Matrix for a 
% Given Problem
%
% Summary: The function calculates the RBF kernel matrix for a given set of 
% range and depth values. The function takes the input vectors `x` and `z` 
% representing the range and depth values respectively. The options struct 
% contains parameters with default values if non-existent or blank fields 
% are encountered.
% 
% SYNTAX:
%   [A,nr] = generateKernel(x, z, options)
%
% INPUTS:
%   x - Vector of range values in the horizontal direction.
%   z - Vector of depth values in the vertical direction.
%   options - Struct containing parameters with the following supported 
%           fields:
%       tau - Resolution of RBF grid with respect to computational grid 
%               (default: 5).
%       eta - Spread of RBF, a positive integer number, denotes how long
%               the effect on neighboring RBF 
%               (default: 4).
%       nouter - Number of RBF layers outside the computational domain 
%               (default: 2).
%       rtype - RBF type: 'global' or 'compact' 
%               (default: 'compact').
%       ltype - Distance type: 'L1', 'L2' or 'LInf' 
%               (default: 'L2').
%
% OUTPUTS:
%   A - RBF kernel matrix with size `n x m`.
%   nr - Row vector with number of RBFs in Z and X direction respectively, 
%       given by [nz, nx].
%
% Author: Ajinkya Kadu
%         Mathematical Institute, Utrecht University, The Netherlands
%
% Date: September 2016

% Validate the number of input arguments
if nargin < 3
    options = [];
end

% Retrieve options
tau     = GetOptions(options, 'tau', 5);
eta     = GetOptions(options, 'eta', 4);
nouter  = GetOptions(options, 'nouter', 2);
rtype   = GetOptions(options, 'rtype', 'compact');
ltype   = GetOptions(options, 'ltype', 'L2');

no = 2 * nouter - 1;
h  = min([x(2) - x(1); z(2) - z(1)]);
hr = tau * h;
[zz, xx] = ndgrid(z, x);

zr = (min(z) - no * hr / 2) : hr : (max(z) + no * hr / 2);
xr = (min(x) - no * hr / 2) : hr : (max(x) + no * hr / 2);
nr = [length(zr), length(xr)]; % Number of RBFs in Z and X direction

[Zc, Xc] = ndgrid(zr, xr);
Zc = Zc(:);
Xc = Xc(:);

nrbf   = length(Zc);
ngridp = length(zz(:));

switch rtype
    case 'global' % Global RBF (Gaussian)
        kernelM = @(r) exp(-r.^2);
        ki = 3.3;
    case 'compact' % Compactly-supported RBF (Wendland C4)
        kernelM = @(r) max(1 - r, 0).^8 .* (32 * r.^3 + 25 * r.^2 + 8 * r + 1);
        ki = 1;
end

% Distance function for RBF
switch ltype
    case 'L1'
        rx = @(x, z) abs(xx(:) - x) + abs(zz(:) - z);
    case 'L2'
        rx = @(x, z) sqrt((xx(:) - x).^2 + (zz(:) - z).^2);
    case 'LInf'
        rx = @(x, z) max(abs(xx(:) - x), abs(zz(:) - z));
end

% Beta values of RBF / Scaling factor of RBF (spread of RBF)
beta = ki * (1 / eta) * (1 / hr);

% RBF kernel matrix
A = zeros(ngridp, nrbf);

for i = 1:nrbf
    A(:, i) = kernelM(beta * rx(Xc(i), Zc(i)));
end

A = sparse(A);

end

