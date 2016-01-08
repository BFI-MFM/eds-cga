function mcompare(s1,s2)
% MCOMPARE :    MCOMPARE ( [ 'file1' ; 'file2' ] , [ 'var1' ; 'var2' ...] )     
%               This optional command plots the relative differences between
%               two different simulations for a list of variables. One plot 
%               is drawn for each variable. The trajectories must have been
%               previously saved by the instruction DYNASAVE. The simulation
%               in file1 serves as the base simulation and the ploted quantity
%               is equal to the difference between the two simulation reported
%               to the first one. If, for a given variable, zero is one of the
%               value of the base simulation, the absolute difference is ploted
%               instead of the relative one.

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

global options_
global nvx nvy x y lag1

ftest(s1,s2) ;

ix = [1-lag1(1):size(x,2)-lag1(1)]' ;
i = [lag1(1):size(ix,1)-lag1(2)+1]' ;

if size(options_.smpl,1) == 1
    error(['DSAMPLE not specified.']) ;
end

if options_.smpl(3) > 0
    if options_.smpl(3) == 2
        if options_.smpl(1)<0 | options_.smpl(2)>size(x,2)-lag1(2)
            error ('Wrong sample.') ;
        end
        i = [options_.smpl(1)+lag1(1):options_.smpl(2)+lag1(1)]' ;
    elseif options_.smpl(3) == 1
        if options_.smpl(1)>size(x,2)-lag1(2)
            error ('Wrong sample.') ;
        end
        i = [lag1(1):options_.smpl(1)+lag1(1)]' ;
    end
end

for k = 1:size(x,1)
    figure ;
    x1 = x(k,i) ;
    y1 = y(k,i) ;
    if nnz(x1) < length(x1)
        plot(ix(i),(y1-x1)) ;
    else
        plot(ix(i),(y1-x1)./x1) ;
    end
    xlabel(['Periods']) ;
    title(['Variable ' s2(k)]) ;
end

return ;
