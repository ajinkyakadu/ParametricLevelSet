function [f,g] = phi(x,d)
%
% Simple inverse problem where forward modeling operator F is Identity matrix
%
% input:
%   x - model vector
%   d - data vector (actual data)
%
% Outputs:
%   f is the misfit value from fitting the computed data from F to actual data 
%   g is the gradient vector
%
% Author: Ajinkya Kadu
%         Mathematical Institute,
%         Utrecht University, The Netherlands
%
% Date : September 2016

F = speye(length(x));
f = 0.5*norm(F*x - d)^2;
g = F'*(F*x - d);

end
 