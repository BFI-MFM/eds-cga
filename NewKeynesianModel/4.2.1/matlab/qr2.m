function [Q,R] = qr2(X)
% This routine performs a qr decomposition of matrix X such that the
% diagonal scalars of the upper-triangular matrix R are positive. If X 
% is a full (column) rank matrix, then R is also the cholesky
% factorization of X'X. This property is needed for the Del Negro 
% & Schorfheides's identification scheme.
% 
% INPUTS 
%   See matlab's documentation for QR decomposition.
%     
% OUTPUTS 
%   See matlab's documentation for QR decomposition.
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

[Q,R] = qr(X);
indx = find(diag(R)<0);
if ~isempty(indx)
    Q(:,indx) = -Q(:,indx);
    R(indx,:) = -R(indx,:);
end