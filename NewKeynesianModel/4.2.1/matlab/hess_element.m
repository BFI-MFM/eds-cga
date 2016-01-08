function d=hess_element(func,element1,element2,args)
% function d=hess_element(func,element1,element2,args)
% returns an entry of the finite differences approximation to the hessian of func
%
% INPUTS
%    func       [function name]    string with name of the function
%    element1   [int]              the indices showing the element within the hessian that should be returned
%    element2   [int]
%    args       [cell array]       arguments provided to func
%
% OUTPUTS
%    d          [double]           the (element1,element2) entry of the hessian
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2010-2011 Dynare Team
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

func = str2func(func);

h=10e-6;
p10 = args;
p01 = args;
m10 = args;
m01 = args;
p11 = args;
m11 = args;
for i=1:size(args,2)
    if i==element1
        p10{i} = p10{i} + h;
        m10{i} = m10{i} - h;

        p11{i} = p11{i} + h;
        m11{i} = m11{i} - h;
    end

    if i==element2
        p01{i} = p01{i} + h;
        m01{i} = m01{i} - h;

        p11{i} = p11{i} + h;
        m11{i} = m11{i} - h;
    end
end

% From Abramowitz and Stegun. Handbook of Mathematical Functions (1965)
% formulas 25.3.24 and 25.3.27 p. 884
if element1==element2
    d = (16*func(p10{:})...
         +16*func(m10{:})...
         -30*func(args{:})...
         -func(p11{:})...
         -func(m11{:}))/(12*h^2);
else
    d = (func(p10{:})...
         +func(m10{:})...
         +func(p01{:})...
         +func(m01{:})...
         -2*func(args{:})...
         -func(p11{:})...
         -func(m11{:}))/(-2*h^2);
end
