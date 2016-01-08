function [g,flag] = moment_function(xparams,sample_moments,dataset,options,parallel)
% Evaluates the moment function of the Simulated Moments Method (discrepancy between sample and
% ).
%
% INPUTS:
%  xparams          [double]  p*1 vector of estimated parameters. 
%  sample_moments   [double]  n*1 vector of sample moments (n>=p).
%  options          [      ]  Structure defining options for SMM.
%  parallel         [      ]  Structure defining the parallel mode settings (optional).
%
% OUTPUTS: 
%  g                [double]  n*1 vector, the gap between simulated and sample moments.
%  flag             [intger]  empty matrix.
%
% SPECIAL REQUIREMENTS
%  The user has to provide a file where the moment conditions are defined.

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

global M_ options_
persistent mainStream mainState
persistent priorObjectiveValue

flag = 1;

if nargin<5
    if isempty(mainStream)
        mainStream = RandStream.getDefaultStream;
        mainState  = mainStream.State;
    else
        mainStream.State = mainState;
    end
end

penalty = 0;
for i=1:options.estimated_parameters.nb
    if ~isnan(options.estimated_parameters.upper_bound(i)) && xparams(i)>options.estimated_parameters.upper_bound(i)
        penalty = penalty + (xparams(i)-options.estimated_parameters.upper_bound(i))^2;
    end
    if ~isnan(options.estimated_parameters.lower_bound(i)) && xparams(i)<options.estimated_parameters.lower_bound(i)
        penalty = penalty + (xparams(i)-options.estimated_parameters.lower_bound(i))^2;
    end
end

if penalty>0
    flag = 0;
    return;
end

save('estimated_parameters.mat','xparams');

% Check for local determinacy of the deterministic steady state.
noprint = options_.noprint; options_.noprint = 1;
[local_determinacy_and_stability,info] = check; options_.noprint = noprint;
if ~local_determinacy_and_stability
    flag = 0;
    return
end

simulated_moments = zeros(size(sample_moments));

% Just to be sure that things don't mess up with persistent variables...
clear perfect_foresight_simulation;

if nargin<5
    for s = 1:options.number_of_simulated_sample
        time_series = extended_path([],options.simulated_sample_size,1);
        data = time_series(dataset.observed_variables_idx,options.burn_in_periods+1:options.simulated_sample_size);
        eval(['tmp = ' options.moments_file_name '(data);'])
        simulated_moments = simulated_moments + tmp;
        simulated_moments = simulated_moments / options.number_of_simulated_sample;
    end
else% parallel mode.
    if ~isunix
        error('The parallel version of SMM estimation is not implemented for non unix platforms!')
    end
    job_number = 1;% Remark. First job is for the master.
    [Junk,hostname] = unix('hostname --fqdn');
    hostname = deblank(hostname);
    for i=1:length(parallel)
        machine = deblank(parallel(i).machine);
        if ~strcmpi(hostname,machine)
            % For the slaves on a remote computer.
            unix(['scp estimated_parameters.mat ' , parallel(i).login , '@' , machine , ':' parallel(i).folder ' > /dev/null']);
        else
            if ~strcmpi(pwd,parallel(i).folder)
                % For the slaves on this computer but not in the same directory as the master.
                unix(['cp estimated_parameters.mat ' , parallel(i).folder]);
            end
        end
        for j=1:parallel(i).number_of_jobs
            if (strcmpi(hostname,machine) && j>1) || ~strcmpi(hostname,machine)  
                job_number = job_number + 1;
                unix(['ssh -A ' parallel(i).login '@' machine ' ./call_matlab_session.sh job' int2str(job_number) '.m &']);
            end
        end
    end
    % Finally the Master do its job
    tStartMasterJob = clock;
    eval('job1;')
    tElapsedMasterJob = etime(clock, tStartMasterJob);
    TimeLimit = tElapsedMasterJob*1.2;
    % Master waits for the  slaves' output... 
    tStart = clock;
    tElapsed = 0;
    while tElapsed<TimeLimit
        if ( length(dir('./intermediary_results_from_master_and_slaves/simulated_moments_slave_*.dat'))==job_number )
            break
        end
        tElapsed = etime(clock, tStart);
    end
    try
        tmp = zeros(length(sample_moments),1);
        for i=1:job_number
            simulated_moments = load(['./intermediary_results_from_master_and_slaves/simulated_moments_slave_' int2str(i) '.dat'],'-ascii');
            tmp = tmp + simulated_moments;
        end
        simulated_moments = tmp / job_number;        
    catch
        flag = 0;
        return
    end
end

g = simulated_moments-sample_moments;