function [x,status]=bicgstab_(func,b,x,tole,kmax,varargin)

% Copyright (C) 2007-2009 Dynare Team
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

status = 0;
r=b-feval(func,x,varargin{:});
rh_0 = r;
rh = r;
rho_0 = 1;
alpha = 1;
w = 1;
v = 0;
p = 0;
k = 0;
rho_1 = rh_0'*r;
tolr = tole*norm(b);

while norm(r) > tolr & k < kmax
    k = k+1;
    beta = (rho_1/rho_0)*(alpha/w);
    p = r+beta*(p-w*v);
    v = feval(func,p,varargin{:});
    alpha = rho_1/(rh_0'*v);
    r = r-alpha*v;
    t = feval(func,r,varargin{:});
    w = (t'*r)/(t'*t);
    rho_0 = rho_1;
    rho_1 = -w*(rh_0'*t);
    x = x+alpha*p+w*r;
    r = r-w*t;
end
if k == kmax
    status = 1;
    warning(sprintf('BICSTABN didn''t converge after %d iterations: norm(r) = %g',kmax,norm(r)));
end