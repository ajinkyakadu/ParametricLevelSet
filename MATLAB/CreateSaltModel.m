function [v] = CreateSaltModel(x,z,v0,options)
%GenerateSaltModel - Generates a velocity model that incorporates a salt 
% deposit embedded within a background model
%
% Summary: The following MATLAB code implements a function that generates 
% a random salt model embedded into a background model. The function, 
% createSaltModel, takes four input arguments and returns a single output 
% argument. The inputs are a range vector x, a depth vector z, a 
% background velocity vector v0, and an options structure, options. The 
% output is a final velocity model containing the salt and the background, 
% represented as a vector v of size n = prod(nx, nz). The options 
% structure allows the user to specify various parameters related to the 
% salt model, including the velocity of the salt, the width and offset of 
% the salt in both the range and depth directions, the number of random 
% points used to form the boundary of the salt, and the random seed used 
% to generate these points. The function employs the inpolygon function 
% to determine the points within the boundary of the salt, which are 
% then assigned the salt velocity.
%
% SYNTAX:
%   [v] = GenerateSaltModel(x,z,v0,options)
%
% INPUTS:
%   x - A vector that represents the range coordinate (horizontal direction)
%   z - A vector that represents the depth coordinate (vertical direction)
%   v0 - Background velocity of the model, represented as a vector with 
%           length n = prod(nx, nz)
%   options - An optional input structure that specifies various 
%           configuration parameters 
%
% OUTPUTS:
%   v - A vector representing the final velocity model, containing both 
%           salt and background, with length n = prod(nx, nz)
%
% Options:
%   options.v1 - Salt velocity in meters per second 
%           (default : 4500)
%   options.xwidth - Percentage width of the salt deposit in the range (x) 
%           direction, expressed as a value between 0 and 1, where 1 
%           indicates complete range coverage 
%           (default : 0.5)
%   options.zwidth - Percentage width of the salt deposit in the depth (z) 
%           direction, expressed as a value between 0 and 1, where 1 
%           indicates complete depth coverage 
%           (default : 0.5)
%   options.xoffset - Offset of the salt deposit in the range (x) 
%           direction, in meters, expressed as a positive or negative value, 
%           with maximum magnitude equal to half of the range 
%           (default : 0)
%   options.zoffset - Offset of the salt deposit in the depth (z) 
%           direction, in meters, expressed as a positive or negative 
%           value, with maximum magnitude equal to half of the depth 
%           (default : 0)
%   options.nrand - The number of random points in the computational 
%           domain used to generate the boundary of the salt deposit 
%           (default : 20)
%   options.randseed - The random seed for generating the random points 
%           that form the boundary of the salt deposit (default : 0)
%
% Author:
%   Ajinkya Kadu
%   Mathematical Institute, Utrecht University, The Netherlands
%
% Date:
%   September 2016

% Validate the number of input arguments
if nargin < 4
    options = [];
end

% Retrieve options
v1      = GetOptions(options, 'v1', 4500);
xwidth  = GetOptions(options, 'xwidth', 0.5);
zwidth  = GetOptions(options, 'zwidth', 0.5);
xoffset = GetOptions(options, 'xoffset', 0);
zoffset = GetOptions(options, 'zoffset', 0);
nrand   = GetOptions(options, 'nrand', 10);
randSeed= GetOptions(options, 'randseed', 0);

% Generate 2D grid of depth and range
[zz,xx] = ndgrid(z,x);

% Initialize velocity model with background velocity
v = v0;

% Set the random number generator seed
rng(randSeed,'twister');

% random points spread over the computational domain
xt = (min(x)+max(x))/2 + xwidth*(max(x)-min(x))*(rand(nrand,1)-0.5) + xoffset;
zt = (min(z)+max(z))/2 + zwidth*(max(z)-min(z))*(rand(nrand,1)-0.5) + zoffset;

% get points on the boundary
k  = boundary(xt,zt);
xv = xt(k);
zv = zt(k);

% Identify points inside the boundary and assign them a salt velocity value
[in,~] = inpolygon(xx,zz,xv,zv);
v(in) = v1;

% Reshape the velocity vector into a column vector
v = v(:);

end

