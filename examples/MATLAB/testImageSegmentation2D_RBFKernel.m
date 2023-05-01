% MATLAB Script for Testing the Parametric Level-Set Method on Simple Inverse Problem
%
% Purpose:
% This script tests the parametric level-set method on a simple inverse 
% problem. The model is generated using 'createSaltModel' which embeds 
% salts into a linearly varying background or a specified background.
% The computational grid used has a size of 10000 m (in range/x-direction) 
% by 3000 m (in depth/z-direction) with a grid spacing of 50 m in each 
% direction. The background has a velocity varying from 1500 m/s to 4000 
% m/s and the salt has a velocity of 4500 m/s.
% The forward modeling operator chosen is an identity matrix.
%
% Author:
% Ajinkya Kadu
% Mathematical Institute, Utrecht University, The Netherlands
%
% Date:
% September 2016


% Clear the command window, clear variables, and close all figures
clc; clearvars; close all;

% Set up the initial grid with dimensions x and z
x = 0:50:10000;
z = 0:50:3000;
[o, d, n] = Grid2Odn(z, x);
[zz, xx] = ngrid(z, x);

% Define the background velocity as a linearly varying value
vBMin = 1500;
vBMax = 4000;
v0 = linspace(vBMin, vBMax, length(z))' * ones(1, length(x));
v0 = v0(:);

% Define the salt velocity as a constant value
v1 = 4500;

SMOptions = [];
SMOptions.v1       = v1;    % Salt velocity
SMOptions.xWidth   = 0.7;   % Width of salt in x direction
SMOptions.zWidth   = 0.6;   % Width of salt in z direction
SMOptions.xOffset  = 0;     % Offset of salt in range (x direction)
SMOptions.zOffset  = 0;     % Offset of salt in depth (z direction)
SMOptions.nRand    = 20;    % Number of points for generating boundary
SMOptions.randSeed = 0;     % Random seed number for generating random points for boundary

% Create the salt model using the defined options
v = CreateSaltModel(x, z, v0, SMOptions);

% Plot the true model
figure(1);
imagesc(x, z, reshape(v, n));
colorbar;
axis equal tight;
title('True Model');

% Generate RBF kernel
KernelOptions.tau   = 5;       % Coarseness of RBF grid relative to computational grid
KernelOptions.eta   = 4;        % Parameter to control RBF spread
KernelOptions.nouter= 2;        % Number of RBF layers outside computational domain
KernelOptions.rtype = 'compact'; % Type of RBF
KernelOptions.ltype = 'L2';     % Distance norms for RBF

[KernelMatrix, NumberOfRBF] = GenerateKernel(x, z, KernelOptions);

% Objective function
ForwardModel = @(x) ComputeMisfitAndGradient(x, v);

PLSOptions.m1     = v1;
PLSOptions.m0     = v0;
PLSOptions.kappa  = 0.1;
ObjectiveFunction = @(x) PLS(x, ForwardModel, KernelMatrix, PLSOptions);

% Initial estimate
InitialEstimate = -1 * ones(NumberOfRBF);
InitialEstimate(end/2 - 2 : end/2 + 2, end/2 - 2 : end/2 + 2) = 1;
InitialEstimate = InitialEstimate(:);

[~, ~, InitialModel] = PLS(InitialEstimate, ForwardModel, KernelMatrix, PLSOptions);

% plot initial estimate
figure(2);
imagesc(x, z, reshape(InitialModel, n));
colorbar; axis equal tight;
title('Initial Model');

% inversion
ModelFittingOptions.maxIter = 100;
ModelFittingOptions.optTol  = 1e-9;
ModelFittingOptions.M       = 5;
FinalEstimate = QGNewton(ObjectiveFunction, InitialEstimate, ModelFittingOptions);

% get final model
PLSOptions.kappa = 0;
[~, ~, FinalModel] = PLS(FinalEstimate, ForwardModel, KernelMatrix, PLSOptions);

% plot final model
figure(3);
imagesc(x, z, reshape(FinalModel, n));
colorbar; axis equal tight;
title('Reconstructed Model');