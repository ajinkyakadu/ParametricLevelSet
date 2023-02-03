% Test script
%
% Testing Parametric Level-Set method on FWI inverse problem
% The model is generated using 'createSaltModel', which embeds salts into a
% linearly varying background (or background of your choice)
% We consider the computational grid of 10000 m (in range/x-direction) by
% 3000 m (in depth/z-direction) with gridspacing of 50 m in each direction.
% The background has velocity varying from 1500 m/s to 4000 m/s. The salt
% has velocity 4500 m/s. The Forward modeling operator is chosen from
% 'SimpleFWI' github repository. Please add it before running
%
% SimpleFWI: https://github.com/TristanvanLeeuwen/SimpleFWI
%
% Author: Ajinkya Kadu
%
% Date : October 2021

clc; clearvars; close all;

% initial setup
dx = 50;
x = 0:dx:10000;
z = 0:dx:3000;
[o,d,n] = grid2odn(z,x);
[zz,xx] = ngrid(z,x);


% background
vbmin = 1500;
vbmax = 4000;
v0 = linspace(vbmin,vbmax,length(z))'*ones(1,length(x));
v0 = v0(:);


% salt
v1 = 4500;
SMoptions.v1        =   v1;             % salt velocity
SMoptions.xwidth    =   0.7;            % width of salt (x-direction)
SMoptions.zwidth    =   0.6;            % width of salt (x-direction) 
SMoptions.xoffset   =   0;              % offset of salt in range (x - direction)
SMoptions.zoffset   =   0;              % offset of salt in depth (z - direction)
SMoptions.nrand     =   20;             % number of points to generate boundary
SMoptions.randseed  =   0;              % random seed number for generating random points for boundary

v = createSaltModel(x,z,v0,SMoptions);
mtrue = 1e6./v.^2;
figure(1); imagesc(x,z,reshape(v,n));colorbar;
axis equal tight;title('true model');


% FWI parameters
dx = 50;
model.f  = 2;                           % frequency in Hz
model.n  = n;                           % model size
model.h  = dx*[1 1];                    % gridspacing
model.xr = 100:2*dx:10000;              % receiver locations (in x dn)
model.zr = 2*dx*ones(1,length(model.xr));% receiver locations (in z dn)
model.xs = 100:4*dx:10000;              % source locations (in x dn)
model.zs = 2*dx*ones(1,length(model.xs));% receiver locations (in z dn)


% measurements (from SimpleFWI github code)
D = F(mtrue(:),model);  


% RBF kernel
Koptions.tau    = 5;        % how coarse the RBF grid should be wrt computational grid
Koptions.eta    = 4;        % parameter to control the spread of RBF
Koptions.nouter = 2;        % RBF layers outside compuational domain 
Koptions.rtype  = 'compact';% RBF type
Koptions.ltype  = 'L2';     % distance norms for RBF

[A,nr] = generateKernel(x,z,Koptions);   


%%% inversion
alpha = 0;                          % regularization parameter
fun   = @(x) misfit(x,D,alpha,model);% FWI misfit function (SimpleFWI code)
PLSoptions.m0    = 1e6./v0.^2;
PLSoptions.m1    = 1e6/v1^2;
PLSoptions.kappa = 0.1;
fh = @(x) PLS(x,fun,A,PLSoptions);

% initial model
x0 = -1*ones(nr);
x0(end/2-2:end/2,end/2-2:end/2) = 1;
x0 = x0(:);
[~,~,m0] = PLS(x0,fun,A,PLSoptions);
figure(2); imagesc(x,z,reshape(sqrt(1e6./m0),n));colorbar;
axis equal tight;title('initial model');

% optimization
mFoptions.maxIter = 100;
mFoptions.optTol = 1e-9;
xf = QGNewton(fh,x0,mFoptions);

% final model
PLSoptions.kappa = 0;
[~,~,mf] = PLS(xf,fun,A,PLSoptions);
vf = sqrt(1e6./mf);
figure(3); imagesc(x,z,reshape(vf,n));colorbar;
axis equal tight;title('reconstructed model');


