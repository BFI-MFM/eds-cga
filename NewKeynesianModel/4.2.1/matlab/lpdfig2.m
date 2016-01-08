function ldens = lpdfig2(x,s,nu)
% Evaluates the logged INVERSE-GAMMA-2 PDF at x.
%
% X ~ IG2(s,nu) if X = inv(Z) where Z ~ G(nu/2,2/s) (Gamma distribution) 
%
% See L. Bauwens, M. Lubrano and J-F. Richard [1999, appendix A] for more details.
%
%
% INPUTS     
%    x     [double]  m*n matrix of locations,
%    s     [double]  m*n matrix or scalar, First INVERSE-GAMMA-2 distribution parameters, 
%    nu    [double]  m*n matrix or scalar, Second INVERSE-GAMMA-2 distribution parameters. 
%
% OUTPUTS
%    ldens [double]  m*n matrix of logged INVERSE-GAMMA-2 densities evaluated at x.
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2004-2009 Dynare Team
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

ldens = -Inf( size(x) ) ;
idx = find( x>0 ) ;    

if length(s)==1
    ldens(idx) = -gammaln(.5*nu) - (.5*nu)*(log(2)-log(s)) - .5*(nu+2)*log(x(idx)) -.5*s./x(idx);
else
    ldens(idx) = -gammaln(.5*nu(idx)) - (.5*nu(idx)).*(log(2)-log(s(idx))) - .5*(nu(idx)+2).*log(x(idx)) -.5*s(idx)./x(idx);
end