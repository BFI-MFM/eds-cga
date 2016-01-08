function x=sylvester3a(x0,a,b,c,d)
% solves iteratively ax+bxc=d

% Copyright (C) 2005-2009 Dynare Team
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

a_1 = inv(a);
b = a_1*b;
d = a_1*d;
e = 1;
iter = 1;
while e > 1e-8 & iter < 500
    x = d-b*x0*c;
    e = max(max(abs(x-x0)));
    x0 = x;
    iter = iter + 1;
end
if iter == 500
    warning('sylvester3a : Only accuracy of %10.8f is achieved after 500 iterations') 
end