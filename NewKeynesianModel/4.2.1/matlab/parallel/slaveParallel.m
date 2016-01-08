function slaveParallel(whoiam,ThisMatlab)
% PARALLEL CONTEXT
% In a parallelization context, this function is launched on slave
% machines, to initialize MATLAB and DYNARE environment and waits for
% instructions sent by the Master. 
% This function is invoked by masterParallel only when the strategy (1),
% i.e. always open, is actived.
%
%
% INPUTS
%  o whoiam [int]       index number of this CPU among all CPUs in the
%                       cluster.
%  o ThisMatlab [int]   index number of this slave machine in the cluster.
%
% OUTPUTS 
%   None  

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

warning off;
diary off;

delete( ['slaveParallel_',int2str(whoiam),'.log']);
diary( ['slaveParallel_',int2str(whoiam),'.log']);


% Configure dynare environment
dynareroot = dynare_config();

% Load input data.
load( ['slaveParallel_input',int2str(whoiam)]);

%Loads fGlobalVar Parallel.
if exist('fGlobalVar'),
    globalVars = fieldnames(fGlobalVar);
    for j=1:length(globalVars),
        eval(['global ',globalVars{j},';']);
        evalin('base',['global ',globalVars{j},';']);
    end
    struct2local(fGlobalVar);
    clear fGlobalVar
    % create global variables in the base workspace as well
    evalin('base',['load( [''slaveParallel_input',int2str(whoiam),'''],''fGlobalVar'')']) ;
    evalin('base','struct2local(fGlobalVar)');
    evalin('base','clear fGlobalVar');
end

t0=clock;
fslave = dir( ['slaveParallel_input',int2str(whoiam),'.mat']);

while (etime(clock,t0)<1200 && ~isempty(fslave)) || ~isempty(dir(['stayalive',int2str(whoiam),'.txt'])),
    if ~isempty(dir(['stayalive',int2str(whoiam),'.txt'])),
        t0=clock;
        delete(['stayalive',int2str(whoiam),'.txt']);
    end
    % I wait for 20 min or while mater asks to exit (i.e. it cancels fslave file)
    pause(1);
    
    fjob = dir(['slaveJob',int2str(whoiam),'.mat']);
    
    if ~isempty(fjob),
        clear fGlobalVar fInputVar fblck nblck fname
        
        while(1)
            Go=0;
            
            Go=fopen(['slaveJob',int2str(whoiam),'.mat']);
            
            if Go>0    
                fclose(Go);
                pause(1);
                load(['slaveJob',int2str(whoiam),'.mat']);
                break
            else
                % Only for testing, will be remouved!
                
                %                if isunix
                %                  E1=fopen('/home/ivano/Works/Errore-slaveParallel.txt','w+');
                %                  fclose(E1);
                %                else            
                %                  E1=fopen('c:\dynare_calcs\Errore-slaveParallel.txt','w+');
                %                  fclose(E1);
                %                end
                
            end
        end
        
        funcName=fname;  % Update global job name.

        if exist('fGlobalVar') && ~isempty (fGlobalVar)
            globalVars = fieldnames(fGlobalVar);
            for j=1:length(globalVars),
                info_whos = whos(globalVars{j});
                if isempty(info_whos) || ~info_whos.global,
                    eval(['global ',globalVars{j},';']);
                    evalin('base',['global ',globalVars{j},';']);
                end
            end
            struct2local(fGlobalVar);
            evalin('base',['load( [''slaveJob',int2str(whoiam),'''],''fGlobalVar'')']);
            evalin('base','struct2local(fGlobalVar)');
            evalin('base','clear fGlobalVar');
        end
        delete(['slaveJob',int2str(whoiam),'.mat']);
        fInputVar.Parallel = Parallel;
        
        % Launch the routine to be run in parallel.
        try,
            tic,
            fOutputVar = feval(fname, fInputVar ,fblck, nblck, whoiam, ThisMatlab);
            toc,
            if isfield(fOutputVar,'OutputFileName'),
                OutputFileName = fOutputVar.OutputFileName;
            else
                OutputFileName = '';
            end

            if(whoiam)

                % Save the output result.
                save([ fname,'_output_',int2str(whoiam),'.mat'],'fOutputVar' );

                % Inform the master that the job is finished, and transfer the output data
                delete(['P_',fname,'_',int2str(whoiam),'End.txt']);
            end

            disp(['Job ',fname,' on CPU ',int2str(whoiam),' completed.']);
            t0 =clock; % Re-set waiting time of 20 mins
        catch ME
            disp(['Job ',fname,' on CPU ',int2str(whoiam),' crashed.']);
            fOutputVar.error = ME;
            save([ fname,'_output_',int2str(whoiam),'.mat'],'fOutputVar' );
            waitbarString = fOutputVar.error.message;
            if Parallel(ThisMatlab).Local,
                waitbarTitle='Local ';
            else
                waitbarTitle=[Parallel(ThisMatlab).ComputerName];
            end
            fMessageStatus(NaN,whoiam,waitbarString, waitbarTitle, Parallel(ThisMatlab));
            delete(['P_',fname,'_',int2str(whoiam),'End.txt']);
            break
            
        end
    end
    fslave = dir( ['slaveParallel_input',int2str(whoiam),'.mat']); % Check if Master asks to exit
end


disp(['slaveParallel on CPU ',int2str(whoiam),' completed.']);
diary off;

delete(['P_slave_',int2str(whoiam),'End.txt']);


exit;
