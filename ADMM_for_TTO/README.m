% Solving a basic truss topology optimization (TTO) problem formulated as
% linear SDP.
%
% To solve a problem, select it from the GEO directory (e.g., t3x3f.geo)
% and do

>> par = kobum('GEO/t3x3.geo');
>> admm_for_tto(par);

% Copyright (c) 2019 Michal Kocvara, m.kocvara@bham.ac.uk
% Last Modified: 15 Dec 2019