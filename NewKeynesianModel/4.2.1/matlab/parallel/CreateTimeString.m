function [TiSt] = CreateTimeString()
% PARALLEL CONTEXT
% In a parallel context, this is a specialized version of clock() function.
%
% INPUTS
% None
%
%  OUTPUTS
%  o TiSt    []   ...  
%
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



TiSt=[];
T=fix(clock);


S1=num2str(T(1));
S2=num2str(T(2));
S3=num2str(T(3));
S4=num2str(T(4));
S5=num2str(T(5));
S6=num2str(T(6));


TiSt=[S1 '-' S2 '-' S3 '-' S4 'h' S5 'm' S6 's'];
