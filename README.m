% ADMM_for_SDP is a collection of MATLAB functions for the solution of
% various SDP and otehr conic problems by the Alternating Direction Method 
% of Multipliers.
% It implements the algorithm published in
%
% Wen, Zaiwen, Donald Goldfarb, and Wotao Yin. "Alternating direction 
% augmented Lagrangian methods for semidefinite programming." Mathematical 
% Programming Computation 2.3-4 (2010): 203-230.
%
% The collection contains the following subdirectories
%
% ADMM_for_SDPA...reads linear SDP problem in SDPA format and solves it
%
% ADMM_for_STANDARD_QP...solves SDP or DNN relaxation of the 
%    standard quadratic programming problem with random data
%
% ADMM_for_TTO...solves the SDP formulation of the basic truss topology
%    optimization problem
%
% All problems can be solved either by 
%    admm_sdp.m ... solves standard linear SDP problem
%    or by
%    admm_dnn.m ... solves a conic problem on the
%    cone of double nonnegative matrices (DNN)
%
% Copyright (c) 2019 Michal Kocvara, m.kocvara@bham.ac.uk
% Last Modified: 15 Dec 2019