function gcompare(s1,s2)
% GCOMPARE :    GCOMPARE ( [ 'file1' ; 'file2' ] , [ 'var1' ; 'var2' ...] )     
%               This optional command plots the trajectories of a list of
%               variables in two different simulations. One plot is drawn
%               for each variable. The trajectories must have been previously
%               saved by the instruction DYNASAVE. The simulation in file1
%               is refered to as the base simulation.

% Copyright (C) 2001-2010 Dynare Team
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

global options_ M_
global nvx nvy x y lag1

ftest(s1,s2) ;

ix = [1-lag1(1):size(x,2)-lag1(1)]' ;
i = [lag1(1):size(ix,1)-lag1(2)+1]' ;

if options_.smpl == 0
    i = [M_.maximum_lag:size(y,2)]' ;
else
    i = [options_.smpl(1):options_.smpl(2)]' ;
end

for k = 1:size(x,1)
    figure ;
    plot (ix(i),x(k,i),ix(i),y(k,i)) ;
    xlabel (['Periods']) ;
    title (['Variable ' s2(k,:)]) ;
    l = min(i) + 1;
    ll = max(i) - 1 ;
    text (l,x(k,l),s1(1,:)) ;
    text (ll,y(k,ll),s1(2,:)) ;
end

% 06/18/01 MJ corrected treatment of options_.smpl
% 06/24/01 MJ removed color specification



