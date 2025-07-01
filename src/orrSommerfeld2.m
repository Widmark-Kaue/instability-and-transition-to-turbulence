function [V,lambda, L, F, y, D, D4] = orrSommerfeld2(omega, Re, U, ddU, opts)
% SPATIALs STABILITY
% Solve Orr-Sommerfeld eigenvalue problem for a specific
% base flow 

arguments
    omega double
    Re double
    U (:,1) double
    ddU (:, 1) double
    opts.H double = 1
    opts.h double = 0
    opts.useSparse logical = false
    opts.numberOfEigenvalues double = 6
    opts.mode (1,:) char {mustBeMember(opts.mode, ["largestabs", ...
        "smallestabs","largestreal",...
        "smallestreal","bothendsreal"])} = "largestabs"
    opts.maxIter double = 100

end

% Get options
H = opts.H;
h = opts.h;
useSparse = opts.useSparse;
K = opts.numberOfEigenvalues;
mode = opts.mode;

N = length(U) + 1;
[D, D2, D4, y] = dirichletNeumannBCs(N, H, h);
    
% Constants

II = eye(size(D));
ZZ = zeros(size(D));

if useSparse
    II = sparse(II);
    ZZ = sparse(ZZ);
    D = sparse(D);
    D2 = sparse(D2);
    D4 = sparse(D4);
    U = sparse(U);
    ddU = sparse(ddU);
end


F0 = -1i*omega*D2 - 1/Re * D4;
F1 = 1i*diag(U)*D2 - 1i*diag(ddU)*II;
F2 = 1i*omega*II + 2/Re * D2;
F3 = -1i*diag(U)*II;
F4 = -1/Re * II;

% set eigenvalue problem
L = [
    ZZ II ZZ ZZ;
    ZZ ZZ II ZZ;
    ZZ ZZ ZZ II;
   -F0 -F1 -F2 -F3;
    ];

F = [
    II ZZ ZZ ZZ;
    ZZ II ZZ ZZ;
    ZZ ZZ II ZZ;
    ZZ ZZ ZZ F4;
    ];

% solve eigenvalue problem
if useSparse
    [V,lambda, flag]=eigs(L,F, K,mode, 'MaxIterations',opts.maxIter);
    
    fprintf('Convergence: %d\n', flag)
    lambda = diag(lambda);
    return
end

[V,lambda]=eig(L,F);
lambda = diag(lambda);
end