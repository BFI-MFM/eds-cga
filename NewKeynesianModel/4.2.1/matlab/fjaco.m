function fjac = fjaco(f,x,varargin)

% FDJAC Computes two-sided finite difference Jacobian
% USAGE
%   fjac = fdjac(f,x,P1,P2,...)
% INPUTS
%   f         : name of function of form fval = f(x)
%   x         : evaluation point
%   P1,P2,... : additional arguments for f (optional)
% OUTPUT
%   fjac      : finite differnce Jacobian
%
% Copyright (C) 2010 Dynare Team
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

ff=feval(f,x,varargin{:});

tol    = eps.^(1/3);
h = tol.*max(abs(x),1);
xh1=x+h; xh0=x-h;
h=xh1-xh0;
fjac = NaN(length(ff),length(x));
for j=1:length(x);
    xx = x;
    xx(j) = xh1(j); f1=feval(f,xx,varargin{:});
    xx(j) = xh0(j); f0=feval(f,xx,varargin{:});
    fjac(:,j) = (f1-f0)/h(j);
end

feval(f,x,varargin{:});
