function [a,b] = ffill(x,ixc,y)
% function [a,b] = ffill(x,ixc,y)
% Makes the horizontal concatenation if x exists
% and fills the matrix with 0 if x and y are not the same size.
%
% INPUTS
%    x:         matrix
%    ixc:       vector of indices
%    y:         matrix
%        
% OUTPUTS
%    a:         concatenation results
%    b:         vector
%        
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2001-2009 Dynare Team
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

xc = size(x,1) ;

if isempty(y)
    b = [ixc; 0];
    a = [x zeros(size(x,1),1)];
else
    yc = size(y,1) ;
    b = [ixc;yc] ;

    if xc > yc
        a = [x [y;zeros(xc-yc,size(y,2))]] ;
    elseif yc > xc
        a = [[x;zeros(yc-xc,size(x,2))] y] ;
    else
        a = [x y] ;
    end

end

% 2001/09/1 MJ corrected for absent lags