function fMessageStatus(prtfrc, njob, waitbarString, waitbarTitle, Parallel)
% PARALLEL CONTEXT
% In parallel context, this function is launched on slave
% machines, and acts as a message passing device for the master machine.

% INPUTS 
% o prtfrc          [double]     fraction of iteration performed
% o njob            [int]        index number of this CPU among all CPUs in the
%                                cluster
% o waitbarString   [char]       running message string to be displayed in the monitor window on master machine 
% o waitbarTitle    [char]       title to be displayed in the monitor window on master machine
% o Parallel        [struct]     options_.parallel(ThisMatlab), i.e. the parallel settings for this slave machine in the cluster.
%
% OUTPUTS 
% o None 

% Copyright (C) 2006-2011 Dynare Team
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

global funcName

if nargin<5,
    Parallel.Local=1;
end

try
    save(['comp_status_',funcName,int2str(njob),'.mat'],'prtfrc','njob','waitbarString','waitbarTitle');
catch  
end

fslave = dir( ['slaveParallel_input',int2str(njob),'.mat']);
if isempty(fslave),
    error('Master asked to break the job');
end

