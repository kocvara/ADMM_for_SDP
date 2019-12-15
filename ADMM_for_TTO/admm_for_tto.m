function admm_for_tto(par)

%
% ADMM_for_TTO - Solves a TTO problem by the ADMM method as 
% described in Z. Wen, D. Goldfarb, and W. Yin. "Alternating direction 
% augmented Lagrangian methods for semidefinite programming." Mathematical
% Programming Computation 2.3-4 (2010): 203-230.
%
% Input file: parameter file generated by function kobum.m
%
% Copyright (c) 2019 Michal Kocvara, m.kocvara@bham.ac.uk
% Last Modified: 15 Dec 2019

m=par.m; n=par.n; n1=par.n1;BI=par.BI;xy=par.xy;
maska=par.maska;ijk=par.ijk; ff=par.f;


nelem = m;
nnod = n1;

gamma=10;
tbar=5.0;

eps=1e-3;%eps=1e-4;

%ff=1.0*circshift(ff,1);

%%
for i=1:m
    x1=xy(ijk(i,2)/2,1); y1=xy(ijk(i,2)/2,2);
    x2=xy(ijk(i,4)/2,1); y2=xy(ijk(i,4)/2,2);
    len(i,1)=sqrt((x1-x2)^2 + (y1-y2)^2);
end

Ahelp=zeros(n,n);t=ones(nelem,1);
for i=1:m
    Ahelp=Ahelp+len(i)*t(i)*BI(i,:)'*BI(i,:);
end
Astiff=Ahelp(maska,maska);

%stiffness matrix
Astiff=sparse(nnod,nnod);
aaa=[];
lb = 1;
for ie=1:nelem
    hel = len(i)*BI(ie,:)'*BI(ie,:);
    ss=sparse(hel(maska,maska));
    aaa = [aaa -ss(:)];
end

am1 = aaa'*aaa + 2.*eye(nelem);
AAT = am1;

clear ss hel am1

C1 = [gamma ff'; ff sparse(nnod,nnod)];
C2 = sparse(nelem,nelem); %C2 = -0.001.*speye(nelem);
C3 = tbar.*speye(nelem);

%% initialize
t = 0.0.*tbar.*ones(nelem,1);
S1 = speye(nnod+1); S2 = speye(nelem); S3 = S2;
X1 = speye(nnod+1); X2 = speye(nelem); X3 = X2; X1o=X1;X2o=X2;X3o=X3;

mu=1; gaga = 0.99; %gaga=.75;

b = -ones(nelem,1);b = -len;

V1m1 = speye(nnod+1);V1m2=V1m1;V2m1 = speye(nelem);V2m2=V2m1;V3m1 = speye(nelem);V3m2=V3m1; VVnormo=1e10;
S1m1 = speye(nnod+1);S1m2=S1m1;S2m1 = speye(nelem);S2m1=S2m1;S3m1 = speye(nelem);S3m2=S3m1;
X1m1 = speye(nnod+1);X1m2=X1m1;X2m1 = speye(nelem);X2m1=X2m1;X3m1 = speye(nelem);X3m2=X3m1;

%%
fprintf('count   pinf          dinf         dgap           VVnorm        err\n');
 count = 1; told=100.*ones(nelem,1); err=1; rho = (1+sqrt(5))/2 - .1; %rho=1;
mup=0;mum=0;VVnorm=1e10;
while err>eps
    
    for iiii=1:1
        %update t
%         Axb = zeros(nelem,1); ASC = zeros(nelem,1);
        XX1 = X1(2:nnod+1,2:nnod+1);SS1 = S1(2:nnod+1,2:nnod+1);
        %ko1 = aaa'*XX1(:); ko2 = aaa'*SS1(:);
        
        Axb = aaa'*XX1(:) - diag(X2) + diag(X3) - b;
        ASC = aaa'*SS1(:) - (diag(S2)-diag(C2)) + (diag(S3)-diag(C3));
        
        t = -AAT\(mu.*Axb + ASC);
        
        %update S
        V1p = C1 - [0 sparse(1,nnod); sparse(nnod,1) reshape(aaa*t,nnod,nnod)];
        V2p = C2 + spdiags(t,0,nelem,nelem);
        V3p = C3 - spdiags(t,0,nelem,nelem);
        V1 = V1p - mu.*X1;
        V2 = V2p - mu.*X2;
        V3 = V3p - mu.*X3;
        
        [evecctor,evalues] = eig(full(V1));
        peval = max(0,evalues);
        S1 = evecctor*peval*evecctor'; S1=sparse(.5.*(S1+S1'));
        S2 = max(0,V2); S3 = max(0,V3);
        
    end
    
    %update X
    X1p = (1/mu).*(S1 - V1);
    X2p = (1/mu).*(S2 - V2);
    X3p = (1/mu).*(S3 - V3);
    
    X1 = (1-rho).*X1 + rho.*X1p;
    X2 = (1-rho).*X2 + rho.*X2p;
    X3 = (1-rho).*X3 + rho.*X3p;
      
    pinf = norm(Axb); pinfs = pinf/(1+norm(b));
    VS = [V1p(:)-S1(:);V2p(:)-S2(:); V3p(:)-S3(:)];
    dinf = sqrt(full(VS'*VS)); dinfs = dinf/(1+norm(ff)+gamma);
    dgap = abs(b'*t - C1(:)'*X1(:)- C2(:)'*X2(:)- C3(:)'*X3(:));
    dgaps = dgap/(1+abs(b'*t)+abs(C1(:)'*X1(:)+C2(:)'*X2(:)+C3(:)'*X3(:)));
    
    err = max(pinfs,max(dinfs,dgaps));
    count = count + 1;
    
    if mod(count,100)==0, 
        pic(par,t); pause(0.00001);
        VVnorm = norm(V1-V1m1,'fro')+norm(V2-V2m1,'fro')+norm(V3-V3m1,'fro');
        if VVnormo<VVnorm-1e-6;
            koc=1;
        end
        fprintf('%3d   %.8f   %.8f   %.8f   %.8f   %.8f\n',count,pinf,dinf,dgap,VVnorm,err);
    end
    V1m2=V1m1;V1m1=V1;V2m2=V2m1;V2m1=V2;V3m2=V3m1;V3m1=V3; VVnormo=VVnorm;
    X1m2=X1m1;X1m1=X1;X2m2=X2m1;X2m1=X2;X3m2=X3m1;X3m1=X3;
    S1m2=S1m1;S1m1=S1;S2m2=S2m1;S2m1=S2;S3m2=S3m1;S3m1=S3;
end
count
err

KK= C1 - [0 sparse(1,nnod); sparse(nnod,1) reshape(aaa*t,nnod,nnod)];
K= -reshape(aaa*t,nnod,nnod);

pic(par,t);