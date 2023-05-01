function [x] = QGNewton(fh, x0, options)
%QGNewton  A simple L-BFGS method with Wolfe linesearch for optimization.
%
%   [xn, info] = QGNewton(fh, x0, options) minimizes an objective function
%   using the L-BFGS method with a Wolfe linesearch strategy.
%
% INPUTS:
%   fh - A function handle to the misfit function. The misfit function must
%        have the form [f, g] = fh(x) where f is the function value and g
%        is the gradient, both with the same size as the input vector x.
%   x0 - The initial guess for the optimization solution.
%   options - An optional structure containing the following fields:
%     maxIter - The maximum number of iterations [default 10].
%     optTol  - The tolerance on the 2-norm of the gradient [default 1e-6].
%     M       - The history size [default 5].
%     fid     - The file ID for output [default 1].
%     write   - A flag indicating whether to save iterates to disk [default 0].
%
% OUTPUTS:
%   xn - The final estimate of the optimization solution.
%
% Author: Tristan van Leeuwen
%         Mathematical Institute, Utrecht University, The Netherlands
% 
% Date: February 2012

if nargin < 3
    options = [];
end

% Parse the options structure.
M         = GetOptions(options, 'M', 5);
fid       = GetOptions(options, 'fid', 1);
itermax   = GetOptions(options, 'maxIter', 10);
tol       = GetOptions(options, 'optTol', 1e-6);
write     = GetOptions(options, 'write', 0);
init_step = GetOptions(options, 'init_step', 1);

% Initialize variables.
n = length(x0);
converged = 0;
iter = 0;
x = x0;
S = zeros(n, 0);
Y = zeros(n, 0);

% Perform initial evaluation.
[f, g] = fh(x);
nfeval = 1;
fprintf(fid, '# iter, # eval, stepsize, f(x), ||g(x)||_2\n');
fprintf(fid, '%6d, %6d, %1.2e, %1.5e, %1.5e\n', iter, nfeval, 1, f, norm(g));
if write
    dlmwrite(['x_' num2str(iter) '.dat'], x);
end

% Main optimization loop.
while ~converged
    % Compute search direction.
    s = B(-g, S, Y);
    p = -(s' * g) / (g' * g);

    if (p < 0)
        fprintf(fid,'Loss of descent, reset history\n');
        S = zeros(n,0);
        Y = zeros(n,0);
        s = B(-g,S,Y);
    end

    % linesearch
    [ft,gt,lambda,lsiter] = WolfeLineSearch(fh,x,f,g,s,init_step);
    nfeval = nfeval + lsiter;

    % update
    xt = x + lambda*s;

    S = [S (xt - x)];
    Y = [Y (gt - g)];

    if (size(S,2) > M)
        S = S(:,end-M+1:end);
        Y = Y(:,end-M+1:end);
    end

    f = ft;
    g = gt;
    x = xt;

    iter = iter + 1;

    fprintf(fid,'%6d, %6d, %1.2e, %1.5e, %1.5e\n',iter,nfeval,lambda,f,norm(g));
    if write
        dlmwrite(['x_' num2str(iter) '.dat'],x);
    end

    % check convergence
    converged = (iter > itermax) || (norm(g) < tol) || (lambda < tol);

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = B(x, S, Y)
% APPLY_LBFGS_INVERSE_HESSIAN_TO_VECTOR
%
% The following function applies the L-BFGS inverse Hessian to a given 
% vector.
%
% INPUTS:
%   x - Column vector of length n.
%   S - Matrix of history of steps with size n x M.
%   Y - Matrix of history of gradient differences with size n x M.
%
% OUTPUTS:
%   z - Column vector of length n, obtained as the result of the application
%       of the L-BFGS inverse Hessian to x.

M     = size(S, 2);

% Initialize variables for calculation.
alpha = zeros(M, 1);
rho   = zeros(M, 1);

% Calculate the values of rho.
for k = 1:M
    rho(k) = 1 / (Y(:, k)' * S(:, k));
end

q = x;

% Perform the first recursion.
for k = M:-1:1
    alpha(k) = rho(k) * S(:, k)' * q;
    q = q - alpha(k) * Y(:, k);
end

% Apply the initial approximation of the Hessian.
if M > 0
    a = (Y(:, end)' * S(:, end)) / (Y(:, end)' * Y(:, end));
else
    a = 1 / norm(x, 1);
end
z = a * q;

% Perform the second recursion.
for k = 1:M
    beta = rho(k) * (Y(:, k)' * z);
    z = z + (alpha(k) - beta) * S(:, k);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ft,gt,lambda,lsiter] = WolfeLineSearch(FunctionHandle, ...
    InitialGuess, InitialFunctionValue, InitialGradient, ...
    SearchDirection, OptionalInitialStep)
% Implements the Simple Wolfe Line Search algorithm, 
% adapted from the source 
% (http://cs.nyu.edu/overton/mstheses/skajaa/msthesis.pdf, Algorithm 3).

lsiter               = 0;
FirstConstant        = 1e-2;
SecondConstant       = 0.9;
TerminationIndicator = false;
LowerBound           = 0;
UpperBound           = Inf;
InitialLambda        = 0.5;

if OptionalInitialStep
    InitialLambda = 0.1 * InitialLambda * norm(InitialGuess) / ...
                norm(SearchDirection);
end

while ~TerminationIndicator
    if UpperBound < Inf
        CurrentLambda = (UpperBound + LowerBound) / 2;
    else
        CurrentLambda = 2 * InitialLambda;
    end

    if lsiter < 10
        [ft,gt] = FunctionHandle(InitialGuess + ...
            (CurrentLambda * SearchDirection));
        lsiter = lsiter + 1;
    else
        CurrentLambda = 0;
        break;
    end

    if (ft > InitialFunctionValue + FirstConstant * CurrentLambda ...
                * (InitialGradient' * SearchDirection))
        UpperBound = CurrentLambda;
    elseif ((gt' * SearchDirection) < SecondConstant * ...
            (InitialGradient' * SearchDirection))
        LowerBound = CurrentLambda;
    else
        TerminationIndicator = true;
    end
end

lambda = CurrentLambda;

end