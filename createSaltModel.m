function [v] = createSaltModel(x,z,v0,options)
%createSaltModel Generates a random salt model embedded into background
%                model
% Summary:
%
% Input:
%   x is a range vector (horizontal direction)
%   z is a depth vector (vertical direction)
%   v0 is Background velocity of model (a vector of size n = prod(nx,nz) )
%
% Output:
%   v is a final velocity model containing salt and background (a vector 
%           of size n = prod(nx,nz)
%
% Supported Input options:
%   v1 - salt velocity in m/s (default : 4500)
%   xwidth - width percentage of the salt in range/x-direction, allowed values 0 to 1, 1 means complete range (default : 0.5)
%   zwidth - width percentage of the salt in depth/z-direction, allowed values 0 to 1, 1 means complete depth (default : 0.5)
%   xoffset - offset of salt in range/x-direction in m, can take +ve and -ve values, maximum is half of range (default : 0)
%   zoffset - offset of salt in depth/z-direction in m, can take +ve and -ve values, maximum is half of depth (default : 0)
%   nrand - number of points in the valid domain to generate the boundary of salt (default : 20)
%   randseed - random seed for generating radom points to form boundary (default : 0)
%
%
% Author: Ajinkya Kadu
%         Mathematical Institute,
%         Utrecht University, The Netherlands
%
% Date : September 2016

if nargin < 4
    options = [];
end

v1      = getoptions(options,'v1',4500);
xwidth  = getoptions(options,'xwidth',0.5);
zwidth  = getoptions(options,'zwidth',0.5);
xoffset = getoptions(options,'xoffset',0);
zoffset = getoptions(options,'zoffset',0);
nrand   = getoptions(options,'nrand',10);
ri      = getoptions(options,'randseed',0);

[zz,xx] = ndgrid(z,x);

v = v0;

rng(ri,'twister');

% random points spread over the computational domain
xt = (min(x)+max(x))/2 + xwidth*(max(x)-min(x))*(rand(nrand,1)-0.5) + xoffset;
zt = (min(z)+max(z))/2 + zwidth*(max(z)-min(z))*(rand(nrand,1)-0.5) + zoffset;

% get points on the boundary
k = boundary(xt,zt);
xv = xt(k);
zv = zt(k);

% locate points inside the boundary and set them to salt velocity
[in,~] = inpolygon(xx,zz,xv,zv);
v(in) = v1;

v = v(:);


end

