function dynareParallelDelete(fname,pname,Parallel)
% PARALLEL CONTEXT
% In a parallel context, this is a specialized version of delete() function.
%
% INPUTS
%  o fname      []   ...
%  o pname      []   ... 
%  o Parallel   []   ...  
%
%  OUTPUTS
%  None
%
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

if nargin ==0,
    disp('dynareParallelDelete(fname)')
    return
end

if nargin ==1,
    pname='';
else
    pname=[pname,filesep];
end

for indPC=1:length(Parallel),
    if isunix
        system(['ssh ',Parallel(indPC).UserName,'@',Parallel(indPC).ComputerName,' rm -f ',Parallel(indPC).RemoteDirectory,'/',pname,fname]);
    else
        delete(['\\',Parallel(indPC).ComputerName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteDirectory,'\',pname,fname]);
    end
    
end
