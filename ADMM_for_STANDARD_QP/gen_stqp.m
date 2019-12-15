function sdpdata = gen_stqp(n); 
% generating data for a standard QP with random data in format used in
% ADMM_for_SDP
%
% Copyright (c) 2019 Michal Kocvara, m.kocvara@bham.ac.uk
% Last Modified: 15 Dec 2019

Q=randn(n,n); Q = Q+Q'; Q = Q'*Q;
%Q = [1 0 1 1 0;0 1 0 1 1;1 0 1 0 1;1 1 0 1 0;0 1 1 0 1];
%Q = diag(diag(ones(n,n)));
E = ones(n,n);

sdpdata.Nx = 1;
sdpdata.Na = 1;
sdpdata.Ng = 0;
sdpdata.d = [];
sdpdata.c = -1;
sdpdata.NaDims = n;
sdpdata.B = [];
sdpdata.A{1} = -sparse(Q);
sdpdata.A{2} = -sparse(E);
sdpdata.Adep{1} = 1;
sdpdata.name = 'random_QP';

