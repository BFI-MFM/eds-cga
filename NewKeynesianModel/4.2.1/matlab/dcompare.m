function dcompare(s1)

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

global options_ nvx nvy x y lag1

ftest(s1,0) ;

i = [lag1(1):size(x,2)-lag1(2)+1]' ;

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

j = bseastr(nvx,nvy) ;

if stop
    return ;
end

z = mean(mean(abs(x(j,i)-y(j,i)))) ;

disp (['The mean absolute difference between set ' s1(1,:) 'and set ' s1(2,:)]) ;
disp (['is : ' num2str(z)]) ;
return ;


