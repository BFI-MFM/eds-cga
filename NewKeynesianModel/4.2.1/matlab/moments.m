function m = moments(X,order)
% Returns the central sample moment of X specified by the positive integer order.
%
% Note that the cross moments are only computed if order=2, in this case the 
% output is a matrix.
%
% INPUTS
%    X      [double]   T*n matrix, where T is the number of observations and n the number of variables.
%    order  [integer]  scalar.
%
% OUTPUTS
%    m      [double]  n*n matrix or n*1 vector of centered moments.
%
% SPECIAL REQUIREMENTS
%    none

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

switch order
  case 1
    m = transpose(mean(X));
  case 2
    m = cov(X);
  otherwise
    if round(order)-order
        error('The second input argument (order) has to be an integer!')
    end
    [T,n] = size(X);
    c = mean(X);
    m = zeros(n,1);
    for i=1:n
        m(i) = mean((X(:,i)-c(i)).^order);
    end
end