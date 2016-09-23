function [A,nr] = generateKernel(x,z,options)
%
% Generates RBF Kernel matrix for given problem
%
% Summary:
%
% Input:
%   x - is a vector with range values (horizontal direction)
%   z - is a vector with depth values (vertical direction)
%   options - is a struct containing parameters (defaults are used for non-existent or blank fields)
%
% Outputs:
%   A is the RBF kernel matrix with size 'n' x 'm'
%   nr is the row vector with number of RBFs in Z and X direction
%   respectively (given by: [nz,nx]).
%
% Supported Input options
%   tau - resolution of RBF grid with respect to computational grid (default : 5)
%   eta - spread of RBF, a positive integer number, denotes how long the
%   effect on the neighboring RBF (default: 4)
%   nouter - RBF layers outside the computational domain (default: 2)
%   rtype - RBF type [ global | (compact) ]
%   ltype - distance type [ L1 | (L2) | LInf ]
%
% Author: Ajinkya Kadu
%         Mathematical Institute,
%         Utrecht University, The Netherlands
%
% Date : September 2016


if nargin < 3
    options = [];
end

tau     = getoptions(options,'tau',5);
eta     = getoptions(options,'eta',4);
nouter  = getoptions(options,'nouter',2);
rtype   = getoptions(options,'rtype','compact');
ltype   = getoptions(options,'ltype','L2');


no  = 2*nouter - 1;

h   = min([x(2)-x(1);z(2)-z(1)]);
hr  = tau*h;
[zz,xx] = ndgrid(z,x);

zr = (min(z)-no*hr/2):hr:(max(z)+no*hr/2);
xr = (min(x)-no*hr/2):hr:(max(x)+no*hr/2);
nr = [length(zr),length(xr)];       % number of RBF in Z and X

[Zc,Xc] = ndgrid(zr,xr);
Zc = Zc(:); Xc = Xc(:);

nrbf    = length(Zc);
ngridp  = length(zz(:));


% RBF type 
switch rtype
    case 'global'          % Global RBF (Gaussian)
        kernelM = @(r) exp(-r.^2);
        ki = 3.3;
    case 'compact'          % Compactly-supported RBF (Wendland C4)
        kernelM = @(r) max(1-r,0).^8.*(32*r.^3 + 25*r.^2 + 8*r + 1);
        ki = 1;
end

% RBF distance function
switch ltype
    case 'L1'
        rx = @(x,z) abs(xx(:)-x) + abs(zz(:)-z);
    case 'L2'
        rx = @(x,z) sqrt((xx(:)-x).^2 + (zz(:)-z).^2);
    case 'LInf'
        rx = @(x,z) max(abs(xx(:)-x),abs(zz(:)-z));
end

% Beta values of RBF / Scaling factor of RBF (spread of RBF)
beta = ki*(1/eta)*(1/hr);

% RBF kernel matrix
A = sparse(ngridp,nrbf);

for i = 1:nrbf
    A(:,i) = kernelM(beta*rx(Xc(i),Zc(i)));
end


end

