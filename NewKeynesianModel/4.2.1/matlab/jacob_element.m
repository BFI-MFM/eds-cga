function d=jacob_element(func,element,args)
% function d=jacob_element(func,element,args)
% returns an entry of the finite differences approximation to the jacobian of func
%
% INPUTS
%    func       [function name]    string with name of the function
%    element    [int]              the index showing the element within the jacobian that should be returned
%    args       [cell array]       arguments provided to func
%
% OUTPUTS
%    d          [double]           jacobian[element]
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
pargs=args;
margs=args;
% length(args) is used instead of size(args, 2) to avoid to transpose column vectors
for i=1:length(args)
    if i==element
        pargs{i} = pargs{i} + h;
        margs{i} = margs{i} - h;
    end
end
d=(func(pargs{:})...
   -func(margs{:}))/(2*h);
