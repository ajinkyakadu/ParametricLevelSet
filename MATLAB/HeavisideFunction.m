function [Heaviside,Delta] = HeavisideFunction(InputVector,Options)
% HeavisideFunction - A function that calculates the Heaviside
% (approximated) and its derivative Delta on a given input vector
%
% Summary:
% The Heaviside function is an approximation of the step function, where
% the step function is defined as H(x) = 1 for x >= 0 and H(x) = 0 for
% x < 0. The Heaviside function is an approximated version of the step
% function, with a smooth transition between the two states. The derivative
% of the Heaviside function is referred to as Delta.
%
% INPUTS:
% InputVector - A real-valued vector of any size
% Options - A structure containing parameters (defaults are used for
%   non-existent or blank fields)
%
% OUTPUTS:
% Heaviside - A vector of Heaviside values, with the same size as the input
%   vector
% Delta - A vector of Delta (derivative of Heaviside) values, with the same
%   size as the input vector
%
% Supported Input options:
% Type - Type of Heaviside [ global | (compact) ]
% Epsi - Heaviside epsilon (default : 0.1)
% Thr - Threshold for level-set, i.e. x = 0 (default : 0)
%
% Author: Ajinkya Kadu
%         Mathematical Institute,
%         Utrecht University, The Netherlands
%
% Date : September 2016

if nargin<2
    Options = [];
end

Type = GetOptions(Options,'Type','compact');
Epsi = GetOptions(Options,'Epsi',0.1);
Thr  = GetOptions(Options,'Thr',0);

InputVector = InputVector - Thr;

switch Type

    case 'global'
        Heaviside = 0.5*(1 + 2/piatan(piInputVector/Epsi));
        Delta = 1./(Epsi*((InputVector.^2*pi^2)/Epsi^2 + 1));

    case 'compact'
        Heaviside = 0*InputVector;
        Delta = 0*InputVector;
        id = find((InputVector < Epsi) & (InputVector > -Epsi));
        Heaviside(id) = 0.5*(1 + InputVector(id)/Epsi + 1/pi*sin(pi*InputVector(id)/Epsi));
        Heaviside(InputVector >= Epsi) = 1;
        Heaviside(InputVector <=-Epsi) = 0;
        Delta(id) = 0.5*(1/Epsi)*(1 + cos(pi*InputVector(id)/Epsi));
end

end