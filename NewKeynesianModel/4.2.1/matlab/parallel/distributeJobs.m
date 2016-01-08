function [nCPU, totCPU, nBlockPerCPU, totSLAVES] = distributeJobs(Parallel, fBlock, nBlock)
% PARALLEL CONTEXT
% In parallel context this function is used to determine the total number of available CPUs,
% and the number of threads to run on each CPU.
%
% INPUTS
%  o Parallel [struct vector]   copy of options_.parallel
%  o fBlock [int]               index number of the first job (e.g. MC iteration or MH block)
%                               (between 1 and nBlock)
%  o nBlock [int]               index number of the last job.
%
% OUTPUT
%  o nBlockPerCPU [int vector]  for each CPU used, indicates the number of
%                               threads run on that CPU
%  o totCPU [int]               total number of CPU used (can be lower than
%                               the number of CPU declared in "Parallel", if
%                               the number of required threads is lower!)
%  o nCPU                       the number of CPU in user format.

% Copyright (C) 2010 Dynare Team
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

totCPU=0;
for j=1:length(Parallel),
    nCPU(j)=length(Parallel(j).CPUnbr);
    totCPU=totCPU+nCPU(j);
end

nCPU=cumsum(nCPU);
offset0 = fBlock-1;
if (nBlock-offset0)>totCPU,
    diff = mod((nBlock-offset0),totCPU);
    nBlockPerCPU(1:diff) = ceil((nBlock-offset0)/totCPU);
    nBlockPerCPU(diff+1:totCPU) = floor((nBlock-offset0)/totCPU);
    totSLAVES=length(Parallel);
else
    nBlockPerCPU(1:nBlock-offset0)=1;
    totCPU = nBlock-offset0;
    totSLAVES = min(find(cumsum(nCPU)>=totCPU));
end
