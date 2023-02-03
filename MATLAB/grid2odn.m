function [grid_centers, grid_spacing, grid_length] = grid2odn(varargin)
% GRID2ODN computes grid information for writing .odn files.
%
% Summary: % The function grid2odn takes in a set of grids, described as 1D
% arrays in each dimension, and returns the grid information required for 
% writing in a .odn file format. The output of the function consists of the
% grid centers, grid spacings, and grid lengths in each dimension. 
% The input arrays are passed as separate arguments, and the function loops
% over each argument to extract the relevant information. The output 
% vectors are initialized to be empty and are appended with the values 
% calculated for each dimension. If the length of a grid in a particular 
% dimension is equal to 1, a grid spacing of 1 is assigned, otherwise the 
% grid spacing is calculated as the difference between the second and first
% element of the input array.
%
% SYNTAX: 
%   [grid_centers, grid_spacing, grid_length] = grid2odn(varargin)
%
% INPUT: 
%   varargin - vectors describing regular grid in each dimension
%
% OUTPUT:
%   grid_centers  - [x(1)        y(1)        ... ] - centers
%   grid_spacing  - [x(2) - x(1) y(2) - y(1) ... ] - gridspacing
%   grid_length   - [length(x)   length(y)   ... ] - grid length
%
% Author: Tristan van Leeuwen
%         Mathematical Institute
%         Utrecht University, The Netherlands
%         
% Date: September 2016

grid_centers = [];
grid_spacing = [];
grid_length = [];

for k=1:length(varargin)
    x = varargin{k};
    grid_centers = [grid_centers x(1)];
    grid_length = [grid_length length(x)];
    if grid_length(end) > 1
        grid_spacing = [grid_spacing x(2) - x(1)];
    else
        grid_spacing = [grid_spacing 1];
    end
end

end

