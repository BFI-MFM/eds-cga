function [sR,sG] = restricted_steadystate(y,x,indx)

% Copyright (C) 2006-2009 Dynare Team
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

global options_ M_ oo_

inde  = options_.steadystate_partial.sseqn;

ss = oo_.steady_state;

ss(indx) = y;

eval(['[R,G] = ' M_.fname '_static(ss, x, M_.params);']);

sR = R(inde);
sG = G(inde,indx);