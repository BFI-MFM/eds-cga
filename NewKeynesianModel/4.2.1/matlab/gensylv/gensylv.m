function [err, E] = gensylv(fake,A,B,C,D)
%function [err, E] = gensylv(fake,A,B,C,D)
% Solves a Sylvester equation.
%
% INPUTS
%   fake     Unused argument (for compatibility with the mex file)
%   A
%   B
%   C
%   D
%    
% OUTPUTS
%   err      [double] scalar: 1 indicates failure, 0 indicates success
%   E
%    
% ALGORITHM
%   none.
%
% SPECIAL REQUIREMENTS
%   none.  

% Copyright (C) 1996-2010 Dynare Team
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

C  = kron(C,C); 
x0 = sylvester3(A,B,C,D);
E  = sylvester3a(x0,A,B,C,D);
err = 0;