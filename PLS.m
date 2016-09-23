function [f,g,m] = PLS(x,funObj,A,options)
%
% PLS Parametric Level-Set function for any kind of F 
%
% Summary: Generates the Parametric Level-Set function for any kind of
% function F on R^{n}. F: R^n to R, PLS: R^m to R. The matrix A,
% called as kernel matrix is defined as A: R^m to R^n. In parametric level
% set, we represent any model vector defined on R^n by m number of RBF
% which is included in kernel matrix A and rephrase the optimization
% problem as finding x in R^m so that the objective defined by F is
% minimized.
% 
% Inputs:
%   x - is a vector (PLS-Parametric Level-Set) defined on R^{m} (for Parametric Level-Set);
%   funObj - is a function handle for solving inverse problem
%   A - RBF kernel matrix with dimensions [ 'n' x 'm' ]
%   options - is a struct containing parameters (defaults are used for non-existent or blank fields)
%
% Outputs:
%   f is the function value from Parametric level-set at PLS vector x.
%   g is the gradient vector of size 'm' at the PLS vector x
%   exitflag returns an exit condition
%   output returns a structure with other information
%
% Supported Input options
%   m0 - Background model vector of size 'n'.
%   kappa - parameter to control width of Heaviside (default : 0.1) 
%   hopt - is a struct containing parameters for Heaviside function
%
% Author: Ajinkya Kadu
%         Mathematical Institute,
%         Utrecht University, The Netherlands
%
% Date : September 2016

if nargin < 4
    options = [];
end

m0    = getoptions(options,'m0',zeros(size(A,1),1));
m1    = getoptions(options,'m1',1);
kappa = getoptions(options,'kappa',0.1);

Ax = A*x;

options.hopt.epsi = kappa*0.5*(max(Ax) - min(Ax));
[h,d] = heavi(Ax,options.hopt);

m = m0.*(1 - h) + h*m1;

[f,g0]  = funObj(m);
g = A'*(g0.*d.*(m1-m0));

end


function [h,d] = heavi(x,options)
%
% heavi Heaviside function (approximated) operating on vector x
%
% Summary
%
% Input:
%   x   - Real-valued vector of any size
%   options - is a struct containing parameters (defaults are used for non-existent or blank fields)
%
% Output:
%   h   = vector of heaviside values of size of x
%   d   = vector of delta (derivative of heaviside) values of size of x
% 
% Supported Input options:
%   type - Type of Heaviside [ global | (compact) ]
%   epsi - Heaviside epsilon (default : 0.01)
%   thr  - Threshold for level-set, i.e x = 0 (default : 0)


if nargin<2
    options = [];
end

type    = getoptions(options,'type','compact');
epsi   = getoptions(options,'epsi',0.1);
thr     = getoptions(options,'thr',0);

x = x - thr;

switch type
    case 'global'
        h = 0.5*(1 + 2/pi*atan(pi*x/epsi));
        d = 1./(epsi*((x.^2*pi^2)/epsi^2 + 1));
        
    case 'compact'
        h = 0*x;
        d = 0*x;
        id = find((x < epsi) & (x > -epsi));
        h(id) = 0.5*(1 + x(id)/epsi + 1/pi*sin(pi*x(id)/epsi));
        h(x >= epsi) = 1;
        h(x <=-epsi) = 0;
        d(id) = 0.5*(1/epsi)*(1 + cos(pi*x(id)/epsi));
        
end

end