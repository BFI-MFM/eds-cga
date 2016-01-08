function [r, g1] = block_mfs_steadystate(y, b, y_all)
% Wrapper around the *_static.m file, for use with dynare_solve,
% when block_mfs option is given to steady.

% Copyright (C) 2009-2010 Dynare Team
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

global M_ oo_

indx = M_.blocksMFS{b};

y_all(indx) = y;
x = [oo_.exo_steady_state; oo_.exo_det_steady_state];

eval(['[r,g1] = ' M_.fname '_static(b, y_all, x, M_.params);']);
g1 = full(g1);