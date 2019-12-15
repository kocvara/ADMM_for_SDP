% admm_sdpa.m is an implementation of ADMM for linear SDP problems
%
% admm_dnn.m is an implementation of ADMM for conic problems with the cone
% double nonnegative matrices
%
% Both functions require input data in the SDPA format
%
% To solve a problem, select it from the "problems" directory 
% (e.g., theta1.dat-s) and do

>> sdpdata = readsdpa('problems/theta1.dat-s');
>> admm_sdp(sdpdata);

% or

>> admm_dnn(sdpdata);

% The maximal nubmer of iterations is set by default to 10000
% This can be changed by an optional second parameter such as

>> admm_sdp(sdpdata,1000);

% Copyright (c) 2019 Michal Kocvara, m.kocvara@bham.ac.uk
% Last Modified: 15 Dec 2019
