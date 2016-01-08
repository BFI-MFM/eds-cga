function x = bseastr(s1,s2)

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

m = size(s1,1) ;
x = zeros(m,1) ;
s1=upper(deblank(s1));
s2=upper(deblank(s2));

for im = 1:m
    key = s1(im,:) ;
    h = size(s2,1) ;
    l = 1 ;
    while l <= h
        mid = round((h+l)/2) ;
        temp = s2(mid,:) ;
        if ~ strcmp(key,temp)
            for i = 1:min(length(key),length(temp))
                if temp(i) > key(i)
                    h = mid - 1 ;
                    break 
                else
                    l = mid + 1 ;
                    break 
                end
            end
        else
            x(im) = mid ;
            break 
        end
    end
end

