function [mu,sigma,offset] = recursive_moments(m0,s0,data,offset)
% Recursive estimation of order one and two moments (expectation and
% covariance matrix). 
% 
% INPUTS 
%   o m0         [double]    (n*1) vector, the prior expectation.
%   o s0         [double]    (n*n) matrix, the prior covariance matrix.
%   o data       [double]    (T*n) matrix.  
%   o offset     [integer]   scalar, number of observation previously
%                            used to compute m0 and s0.
% OUTPUTS 
%   o mu         [double]    (n*1) vector, the posterior expectation. 
%   o sigma      [double]    (n*n) matrix, the posterior covariance matrix.
%   o offset     [integer]   = offset + T.
%
% ALGORITHM 
%   None.       
%
% SPECIAL REQUIREMENTS
%   None.

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

[T,n] = size(data);

for t = 1:T
    tt = t+offset;
    m1 = m0 + (data(t,:)'-m0)/tt;
    qq = m1*m1';
    s1 = s0 + ( (data(t,:)'*data(t,:)-qq-s0) + (tt-1)*(m0*m0'-qq') )/tt;
    m0 = m1;
    s0 = s1;
end

mu = m1; 
sigma = s1;
offset = offset+T;