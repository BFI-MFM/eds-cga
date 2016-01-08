function [A,B] = transition_matrix(dr, varargin)
% function [A,B] = transition_matrix(dr, varargin)
% Makes transition matrices out of ghx and ghu
%
% INPUTS
%    dr:        structure of decision rules for stochastic simulations
% varargin:     {1}: M_
%
% OUTPUTS
%    A:         matrix of effects of predetermined variables in linear solution (ghx)
%    B:         matrix of effects of shocks in linear solution (ghu)
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2009 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if(length(varargin)<=0)
    global M_
else
    M_=varargin{1};
end;

exo_nbr = M_.exo_nbr;
ykmin_ = M_.maximum_endo_lag;

nx = size(dr.ghx,2);
kstate = dr.kstate;
ikx = [dr.nstatic+1:dr.nstatic+dr.npred];

A = zeros(nx,nx);
k0 = kstate(find(kstate(:,2) <= ykmin_+1),:);
i0 = find(k0(:,2) == ykmin_+1);
A(i0,:) = dr.ghx(ikx,:);
B = zeros(nx,exo_nbr);
if(isfield(dr,'ghu'))
    B(i0,:) = dr.ghu(ikx,:);
end;
for i=ykmin_:-1:2
    i1 = find(k0(:,2) == i);
    n1 = size(i1,1);
    j = zeros(n1,1);
    for j1 = 1:n1
        j(j1) = find(k0(i0,1)==k0(i1(j1),1));
    end
    A(i1,i0(j))=eye(n1);
    i0 = i1;
end
