function admm_sdp(sdpdata,maxiter)
%
% ADMM_GEN - Solves general linear SDP problems by the ADMM method as 
% described in Z. Wen, D. Goldfarb, and W. Yin. "Alternating direction 
% augmented Lagrangian methods for semidefinite programming." Mathematical
% Programming Computation 2.3-4 (2010): 203-230.
%
% Input file: Matlab stucture in Penlab format
%
% Elements of the structure
%   name ... filename of the input file
%   Nx ..... number of primal variables
%   Na ..... number of linear matrix inequalities (or diagonal blocks of the
%            matrix constraint)
%   Ng ..... number of linear inequalitites
%   B ...... matrix defining the linear inequality constraints Bx<=d
%            dimension Ng x Nx
%   d ...... rhs for linear constraints
%   c ...... dim (Nx,1), coefficients of the linear objective function
%   NaDims . vector of sizes of matrix constraints (diagonal blocks)
%   A ...... cell array (matrix) of A{k,l} for k=1,...,Na matrix constraint
%            for l=1 ~ absolute term, l=2..Nx+1 coeficient matrices
%            (some of them might be empty)
%   Adep ... dependency list for each matrix constraint
%
% Use the function "readsdpa" to convert SDPA input file into this format:
%
% sdpdata = readsdpa('theta1.dat-s');
%
% Copyright (c) 2019 Michal Kocvara, m.kocvara@bham.ac.uk
% Last Modified: 15 Dec 2019

if nargin<2
   maxiter = 1e4;
end

%% Read SDP input data
n = sdpdata.Nx;
ncon = sdpdata.Na;
b = sdpdata.c;

for icon=1:ncon
    C{icon} = -sdpdata.A{icon,1};
    m(icon) = size(C{icon},1);
end

AAT = sparse(n,n);
for icon=1:ncon
    aaa{icon} = [];
    for i=1:n
        A{icon,i} = sdpdata.A{icon,i+1};
        aaa{icon} = [aaa{icon} A{icon,i}(:)];
    end
    AAT = AAT + aaa{icon}'*aaa{icon};
end

fprintf('================================= A D M M  for  S D P ===============================\n');
fprintf('Number of LMI Constraints: %3d\n',ncon);
fprintf('Number of Variables: %3d\n',n);
fprintf('Maximal Constraint Size: %3d\n',max(m));
fprintf('Problem Name: "%s"\n',sdpdata.name);

%% initializing variables

y = ones(n,1);
for icon=1:ncon
    S{icon} = speye(m(icon));
    X{icon} = speye(m(icon));
end

%% setting parameters
fprintf('-------------------------------------------------------------------------------------\n');
fprintf(' iter    p-infeas     d-infeas       d-gap           mu         error      objective\n');
fprintf('-------------------------------------------------------------------------------------\n');

eps=1e-3;   eps=1e-6; % final precision required
rho = (1+sqrt(5))/2 - .5; %rho=1.;   % step-length for the multiplier X

mu=500; % Augmented Lagrangian penalty parameter
% parameters for update of mu
gamma = 0.5; mu_min=1e-4; mu_max = 1e4; eta1 = 100; eta2 = 100; h4 = 100;

count = 1;  err=1; it_pinf = 0; it_dinf = 0;

tic
%% ADMM main loop
while err>eps
    
    %update y
    Axb = sparse(n,1); ASC = sparse(n,1);
    for icon=1:ncon
        Axb = Axb + aaa{icon}'*X{icon}(:);
        ASC = ASC + aaa{icon}'*(S{icon}(:)-C{icon}(:));
    end
    y = -AAT\(mu.*(Axb-b) + ASC);
    
    %update S
    for icon=1:ncon
        Vp{icon} = C{icon} - [reshape(aaa{icon}*y,m(icon),m(icon))];
        V{icon} = Vp{icon} - mu.*X{icon};      
        [evecctor,evalues] = eig(full(V{icon}));
        peval = max(0,evalues);
        S{icon} = evecctor*peval*evecctor';
        S{icon}=sparse(.5.*(S{icon}+S{icon}'));
    end
    
    %update X
    for icon=1:ncon
        Xp = (1/mu).*(S{icon} - V{icon});
        X{icon} = (1-rho).*X{icon} + rho.*Xp;
    end
    
    % calculate current error
    dinf = 0; dinfs = 0; dgap = 0; dgaps = 0;
    pinf = norm(Axb-b); pinfs = pinf/(1+norm(b));
    for icon=1:ncon
        VS = [Vp{icon}(:)-S{icon}(:)];
        dinfi = sqrt(full(VS'*VS)); dinfsi = dinfi/(1+norm(C{icon},1));
        dinf = dinf + dinfi; dinfs = dinfs + dinfsi;
        dgapi = C{icon}(:)'*X{icon}(:);
        dgapsi = abs(C{icon}(:)'*X{icon}(:));
        dgap = dgap + dgapi; dgaps = dgaps + dgapsi;
    end
    dgap = abs(b'*y - dgap); dgaps = dgap/(1+abs(b'*y)+dgaps);
    
    err = max(pinfs,max(dinfs,dgaps));
    count = count + 1;
    
    % update penalty parameter mu
    if pinf+dinf>2
        if pinf/dinf < eta1
            it_pinf = it_pinf + 1; it_dinf = 0;
            if it_pinf > h4
                mu = max(gamma*mu,mu_min); it_pinf = 0;
            end
        else
            if pinf/dinf > eta2
                it_dinf = it_dinf + 1; it_finf = 0;
                if it_dinf > h4
                    mu = min(mu/gamma,mu_max); it_dinf = 0;
                end
            end
        end
    end
    
    if mod(count,100)==0,
        fprintf('%5d   %.8f   %.8f   %.8f   %.8f   %.6f   %.8f\n',count,pinf,dinf,dgap,mu,err,b'*y);
    end
    
    if count > maxiter, break, end
end
elapsedtime = toc;

fprintf('%5d   %.8f   %.8f   %.8f   %.8f   %.6f   %.8f\n',count,pinf,dinf,dgap,mu,err,b'*y);
fprintf('-------------------------------------------------------------------------------------\n');
lambda_min = min(eig(X{1})); min_el = min(min(X{1}));
fprintf('Minimal eigenvalue of X: %4.2e; Minimal element of X:  %4.2e\n',lambda_min,min_el);
fprintf('-------------------------------------------------------------------------------------\n');
fprintf('Total ADMM iterations: %3d; Final precision: %4.2e; CPU Time %2.2fsec\n',count,err,elapsedtime);
fprintf('=====================================================================================\n');





