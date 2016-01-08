function d = bksup0(c,ny,jcf,iyf,icf,periods)
% Solves deterministic models recursively by backsubstitution for one lead/lag
%
% INPUTS
%    ny:             number of endogenous variables
%    jcf:            variables index forward
%    
% OUTPUTS
%    d:              vector of backsubstitution results
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

ir  = ((periods-2)*ny+1):(ny+(periods-2)*ny);
irf = iyf+(periods-1)*ny ;

for i = 2:periods
    c(ir,jcf) = c(ir,jcf)-c(ir,icf)*c(irf,jcf);
    ir = ir-ny;
    irf = irf-ny;
end

d = c(:,jcf);