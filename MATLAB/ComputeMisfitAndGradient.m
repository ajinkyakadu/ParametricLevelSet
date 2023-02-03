function [FunctionValue, GradientValue] = ComputeMisfitAndGradient(ModelVector, ObservedData)
%ComputeMisfitAndGradient This function computes the misfit and gradient 
% values for a simple inverse problem where the forward modeling operator 
% is represented by an identity matrix.
%
% INPUTS:
%   ModelVector - A vector representing the model parameters
%   ObservedData - A vector representing the observed data
%
% OUTPUTS:
%   FunctionValue - A scalar value representing the misfit value, which is 
%       the difference between the computed data from the forward modeling 
%       operator and the observed data
%   GradientValue - A vector representing the gradient of the misfit function
%
% Author: Ajinkya Kadu
%         Mathematical Institute,
%         Utrecht University, The Netherlands
%
% Date: September 2016

residual      = ModelVector - ObservedData;
FunctionValue = 0.5 * norm(residual)^2;
GradientValue = residual;

end