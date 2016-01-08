function storeGlobalVars(fname,append)
% PARALLEL CONTEXT
% In a parallel context, this function stores all global vars in structure 
% fGlobalVar and saves it in the file fname.mat
%
% INPUTS
% fname  [str]         name of the file
%
% append []                        flag to append globals to the storage file
%
% OUTPUTS
% None
%
%
% Copyright (C) 2009-2011 Dynare Team
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


GlobalNames = who('global');
% varlist = '';
for j=1:length(GlobalNames);
    eval(['global ',GlobalNames{j},';']);
    eval(['fGlobalVar.',GlobalNames{j},'=',GlobalNames{j},';']);
end

if nargin<2,
    save(fname,'fGlobalVar');
else
    save(fname,'fGlobalVar','-append');
end