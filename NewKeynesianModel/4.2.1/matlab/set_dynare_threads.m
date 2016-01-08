function set_dynare_threads(mexname,n)
% This function sets the number of threads used by some MEX files when compiled
% with OpenMP support (i.e with --enable-openmp is given to configure) or any 
% other parallel library.
%
% INPUTS 
%  o mexname  [string]    Name of the mex file.  
%  o n        [integer]   scalar specifying the number of threads to be used.
%
% OUTPUTS 
%  none.

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
global options_

if ~ischar(mexname)
    error('set_dynare_threads:: First argument has to be a string!')
end

if ~isint(n)
    error('set_dynare_threads:: Second argument has to be an integer!')
end

switch mexname
  case 'A_times_B_kronecker_C'
    options_.threads.kronecker.A_times_B_kronecker_C = n;
  case 'sparse_hessian_times_B_kronecker_C'
    options_.threads.kronecker.sparse_hessian_times_B_kronecker_C = n;
  otherwise
    message = [ mexname ' is not a known parallel mex file.' ];
    message_id  = 'Dynare:Threads:UnknownParallelMex';
    warning(message_id,message);
end