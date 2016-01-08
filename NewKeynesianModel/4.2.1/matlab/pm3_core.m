function myoutput=pm3_core(myinputs,fpar,nvar,whoiam, ThisMatlab)

% PARALLEL CONTEXT
% Core functionality for pm3.m function, which can be parallelized.

% INPUTS 
% See the comment in random_walk_metropolis_hastings_core.m funtion.

% OUTPUTS
% o myoutput  [struc]
%
%
% ALGORITHM 
%   Portion of McMCDiagnostics.m function.       
%
% SPECIAL REQUIREMENTS.
%   None.

% Copyright (C) 2007-2011 Dynare Team
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

if nargin<4,
    whoiam=0;
end

% Reshape 'myinputs' for local computation.
% In order to avoid confusion in the name space, the instruction struct2local(myinputs) is replaced by:

tit1=myinputs.tit1;
nn=myinputs.nn;
n2=myinputs.n2;
Distrib=myinputs.Distrib;
varlist=myinputs.varlist;
MaxNumberOfPlotsPerFigure=myinputs.MaxNumberOfPlotsPerFigure;
name3=myinputs.name3;
tit3=myinputs.tit3;
Mean=myinputs.Mean;

if whoiam
    Parallel=myinputs.Parallel;
end


global options_ M_ oo_


if whoiam
    waitbarString = ['Parallel plots pm3 ...'];
    if Parallel(ThisMatlab).Local,
        waitbarTitle=['Local '];
    else
        waitbarTitle=[Parallel(ThisMatlab).ComputerName];
    end        
    fMessageStatus(0,whoiam,waitbarString, waitbarTitle, Parallel(ThisMatlab));   
end



figunumber = 0;
subplotnum = 0;
hh = figure('Name',[tit1 ' ' int2str(figunumber+1)]);
RemoteFlag = 0;
if whoiam,
    if Parallel(ThisMatlab).Local ==0
        RemoteFlag=1;
    end
end

OutputFileName = {};

for i=fpar:nvar
    if max(abs(Mean(:,i))) > 10^(-6)
        subplotnum = subplotnum+1;
        set(0,'CurrentFigure',hh);
        subplot(nn,nn,subplotnum);
        plot([1 n2],[0 0],'-r','linewidth',0.5);
        hold on
        for k = 1:9
            plot(1:n2,squeeze(Distrib(k,:,i)),'-g','linewidth',0.5);
        end
        plot(1:n2,Mean(:,i),'-k','linewidth',1);
        xlim([1 n2]);
        hold off;
        name = deblank(varlist(i,:));
        title(name,'Interpreter','none')
    end
    
    if whoiam,
        if Parallel(ThisMatlab).Local==0
            DirectoryName = CheckPath('Output');
        end
    end
    
    if subplotnum == MaxNumberOfPlotsPerFigure | i == nvar
        eval(['print -depsc2 ' M_.dname '/Output/'  M_.fname '_' name3 '_' deblank(tit3(i,:)) '.eps' ]);
        if ~exist('OCTAVE_VERSION')
            eval(['print -dpdf ' M_.dname '/Output/' M_.fname  '_' name3 '_' deblank(tit3(i,:))]);
            saveas(hh,[M_.dname '/Output/' M_.fname '_' name3 '_' deblank(tit3(i,:)) '.fig']);
        end
        if RemoteFlag==1,
            OutputFileName = [OutputFileName; {[M_.dname, filesep, 'Output',filesep], [M_.fname '_' name3 '_' deblank(tit3(i,:)) '.*']}];
        end
        if options_.nograph, close(hh), end
        subplotnum = 0;
        figunumber = figunumber+1;
        if (i ~= nvar)
            hh = figure('Name',[name3 ' ' int2str(figunumber+1)]);
        end
    end
    
    if whoiam,
        waitbarString = [ 'Variable ' int2str(i) '/' int2str(nvar) ' done.'];
        fMessageStatus((i-fpar+1)/(nvar-fpar+1),whoiam,waitbarString, waitbarTitle, Parallel(ThisMatlab));
    end
    
    
end


myoutput.OutputFileName=OutputFileName;
