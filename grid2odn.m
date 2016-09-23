function [o,d,n] = grid2odn(varargin)
% produces grid info for writing .odn files.
%
% use:
%   [o,d,n] = grid2odn(x,y,...)
%
% input:
%   {x,y} - vectors describing regular grid in each dimension
%
% output:
%   o - [x(1)        y(1)        ... ] - centers
%   d - [x(2) - x(1) y(2) - y(1) ... ] - gridspacing
%   n - [length(x)   length(y)   ... ] - grid length
%
% Author: Tristan van Leeuwen
%         Mathematical Institute
%         Utrecht University, The Netherlands
%         
% Date: September 2016


o = [];
d = [];
n = [];

for k=1:length(varargin)
    x = varargin{k};
    o = [o x(1)];
    n = [n length(x)];
    if n(end)>1
        d = [d x(2)-x(1)];
    else
        d = [d 1];
    end
end