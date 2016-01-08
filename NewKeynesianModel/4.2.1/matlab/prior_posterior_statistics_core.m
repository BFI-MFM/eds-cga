function myoutput=prior_posterior_statistics_core(myinputs,fpar,B,whoiam, ThisMatlab)
% PARALLEL CONTEXT
% Core functionality for prior_posterior.m function, which can be parallelized.
% See also the comment in random_walk_metropolis_hastings_core.m funtion.
%
% INPUTS
%   See See the comment in random_walk_metropolis_hastings_core.m funtion.
%
% OUTPUTS
% o myoutput  [struc]
%  Contained OutputFileName_smooth;
%                          _update;
%                          _inno;
%                          _error;
%                          _filter_step_ahead;
%                          _param;
%                          _forc_mean;
%                          _forc_point
%
% ALGORITHM
%   Portion of prior_posterior.m function.
% This file is part of Dynare.
%
% SPECIAL REQUIREMENTS.
%   None.

% Copyright (C) 2005-2011 Dynare Team
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

global options_ oo_ M_ bayestopt_ estim_params_

if nargin<4,
    whoiam=0;
end

global options_ oo_ M_ bayestopt_ estim_params_


% Reshape 'myinputs' for local computation.
% In order to avoid confusion in the name space, the instruction struct2local(myinputs) is replaced by:

type=myinputs.type;
run_smoother=myinputs.run_smoother;
gend=myinputs.gend;
Y=myinputs.Y;
data_index=myinputs.data_index;
missing_value=myinputs.missing_value;
varobs=myinputs.varobs;
irun=myinputs.irun;
endo_nbr=myinputs.endo_nbr;
nvn=myinputs.nvn;
naK=myinputs.naK;
horizon=myinputs.horizon;
iendo=myinputs.iendo;
if horizon
    i_last_obs=myinputs.i_last_obs;
    IdObs=myinputs.IdObs;
    MAX_nforc1=myinputs.MAX_nforc1;
    MAX_nforc2=myinputs.MAX_nforc2;
end
if naK
    MAX_naK=myinputs.MAX_naK;
end

exo_nbr=myinputs.exo_nbr;
maxlag=myinputs.maxlag;
MAX_nsmoo=myinputs.MAX_nsmoo;
MAX_ninno=myinputs.MAX_ninno;
MAX_nerro = myinputs.MAX_nerro;
MAX_nruns=myinputs.MAX_nruns;
MAX_momentsno = myinputs.MAX_momentsno;
ifil=myinputs.ifil;

if ~strcmpi(type,'prior'),
    x=myinputs.x;
    if strcmpi(type,'posterior'),
        logpost=myinputs.logpost;
    end
end
if whoiam
    Parallel=myinputs.Parallel;
end

% DirectoryName = myinputs.DirectoryName;
if strcmpi(type,'posterior')
    DirectoryName = CheckPath('metropolis');
elseif strcmpi(type,'gsa')
    if options_.opt_gsa.pprior
        DirectoryName = CheckPath(['GSA',filesep,'prior']);
    else
        DirectoryName = CheckPath(['GSA',filesep,'mc']);
    end
elseif strcmpi(type,'prior')
    DirectoryName = CheckPath('prior');
end

RemoteFlag = 0;
if whoiam
    for j=1:length(Parallel),
        if Parallel(j).Local==0,
            RemoteFlag = 1;
        end
    end
    ifil=ifil(:,whoiam);
    waitbarString = ['Please wait... ',type,' subdraws (' int2str(fpar) 'of' int2str(B) ')...'];
    if Parallel(ThisMatlab).Local,
        waitbarTitle=['Local '];
    else
        waitbarTitle=[Parallel(ThisMatlab).ComputerName];
    end
    fMessageStatus(0,whoiam,waitbarString, waitbarTitle, Parallel(ThisMatlab));
else
    if exist('OCTAVE_VERSION')
        diary off;
    else
        h = waitbar(0,['Taking ',type,' subdraws...']);
    end

end

if RemoteFlag==1,
    OutputFileName_smooth = {};
    OutputFileName_update = {};
    OutputFileName_inno = {};
    OutputFileName_error = {};
    OutputFileName_filter_step_ahead = {};
    OutputFileName_param = {};
    OutputFileName_forc_mean = {};
    OutputFileName_forc_point = {};
    % OutputFileName_moments = {};
end

for b=fpar:B

    %    [deep, logpo] = GetOneDraw(type);
    %    set_all_parameters(deep);
    %    dr = resol(oo_.steady_state,0);
    if strcmpi(type,'prior')

        [deep, logpo] = GetOneDraw(type);

    else
        deep = x(b,:);
        if strcmpi(type,'posterior')
            logpo = logpost(b);
        else
            logpo = evaluate_posterior_kernel(deep');
        end
    end
    set_all_parameters(deep);
    [dr,info] = resol(oo_.steady_state,0);

    if run_smoother
        [alphahat,etahat,epsilonhat,alphatilde,SteadyState,trend_coeff,aK] = ...
            DsgeSmoother(deep,gend,Y,data_index,missing_value);

        if options_.loglinear
            stock_smooth(dr.order_var,:,irun(1)) = alphahat(1:endo_nbr,:)+ ...
                repmat(log(SteadyState(dr.order_var)),1,gend);
            stock_update(dr.order_var,:,irun(1)) = alphatilde(1:endo_nbr,:)+ ...
                repmat(log(SteadyState(dr.order_var)),1,gend);
        else
            stock_smooth(dr.order_var,:,irun(1)) = alphahat(1:endo_nbr,:)+ ...
                repmat(SteadyState(dr.order_var),1,gend);
            stock_update(dr.order_var,:,irun(1)) = alphatilde(1:endo_nbr,:)+ ...
                repmat(SteadyState(dr.order_var),1,gend);
        end
        stock_innov(:,:,irun(2))  = etahat;
        if nvn
            stock_error(:,:,irun(3))  = epsilonhat;
        end
        if naK
            stock_filter_step_ahead(:,dr.order_var,:,irun(4)) = aK(options_.filter_step_ahead,1:endo_nbr,:);
        end

        if horizon
            yyyy = alphahat(iendo,i_last_obs);
            yf = forcst2a(yyyy,dr,zeros(horizon,exo_nbr));
            if options_.prefilter == 1
                yf(:,IdObs) = yf(:,IdObs)+repmat(bayestopt_.mean_varobs', ...
                                                 horizon+maxlag,1);
            end
            yf(:,IdObs) = yf(:,IdObs)+(gend+[1-maxlag:horizon]')*trend_coeff';
            if options_.loglinear == 1
                yf = yf+repmat(log(SteadyState'),horizon+maxlag,1);
            else
                yf = yf+repmat(SteadyState',horizon+maxlag,1);
            end
            yf1 = forcst2(yyyy,horizon,dr,1);
            if options_.prefilter == 1
                yf1(:,IdObs,:) = yf1(:,IdObs,:)+ ...
                    repmat(bayestopt_.mean_varobs',[horizon+maxlag,1,1]);
            end
            yf1(:,IdObs,:) = yf1(:,IdObs,:)+repmat((gend+[1-maxlag:horizon]')* ...
                                                   trend_coeff',[1,1,1]);
            if options_.loglinear == 1
                yf1 = yf1 + repmat(log(SteadyState'),[horizon+maxlag,1,1]);
            else
                yf1 = yf1 + repmat(SteadyState',[horizon+maxlag,1,1]);
            end

            stock_forcst_mean(:,:,irun(6)) = yf';
            stock_forcst_point(:,:,irun(7)) = yf1';
        end

    end
    stock_param(irun(5),:) = deep;
    stock_logpo(irun(5),1) = logpo;
    stock_ys(irun(5),:) = SteadyState';

    irun = irun +  ones(7,1);


    if irun(1) > MAX_nsmoo || b == B
        stock = stock_smooth(:,:,1:irun(1)-1);
        ifil(1) = ifil(1) + 1;
        save([DirectoryName '/' M_.fname '_smooth' int2str(ifil(1)) '.mat'],'stock');

        stock = stock_update(:,:,1:irun(1)-1);
        save([DirectoryName '/' M_.fname '_update' int2str(ifil(1)) '.mat'],'stock');
        if RemoteFlag==1,
            OutputFileName_smooth = [OutputFileName_smooth; {[DirectoryName filesep], [M_.fname '_smooth' int2str(ifil(1)) '.mat']}];
            OutputFileName_update = [OutputFileName_update; {[DirectoryName filesep], [M_.fname '_update' int2str(ifil(1)) '.mat']}];
        end
        irun(1) = 1;
    end

    if irun(2) > MAX_ninno || b == B
        stock = stock_innov(:,:,1:irun(2)-1);
        ifil(2) = ifil(2) + 1;
        save([DirectoryName '/' M_.fname '_inno' int2str(ifil(2)) '.mat'],'stock');
        if RemoteFlag==1,
            OutputFileName_inno = [OutputFileName_inno; {[DirectoryName filesep], [M_.fname '_inno' int2str(ifil(2)) '.mat']}];
        end
        irun(2) = 1;
    end

    if nvn && (irun(3) > MAX_nerro || b == B)
        stock = stock_error(:,:,1:irun(3)-1);
        ifil(3) = ifil(3) + 1;
        save([DirectoryName '/' M_.fname '_error' int2str(ifil(3)) '.mat'],'stock');
        if RemoteFlag==1,
            OutputFileName_error = [OutputFileName_error; {[DirectoryName filesep], [M_.fname '_error' int2str(ifil(3)) '.mat']}];
        end
        irun(3) = 1;
    end

    if naK && (irun(4) > MAX_naK || b == B)
        stock = stock_filter_step_ahead(:,:,:,1:irun(4)-1);
        ifil(4) = ifil(4) + 1;
        save([DirectoryName '/' M_.fname '_filter_step_ahead' int2str(ifil(4)) '.mat'],'stock');
        if RemoteFlag==1,
            OutputFileName_filter_step_ahead = [OutputFileName_filter_step_ahead; {[DirectoryName filesep], [M_.fname '_filter_step_ahead' int2str(ifil(4)) '.mat']}];
        end
        irun(4) = 1;
    end

    if irun(5) > MAX_nruns || b == B
        stock = stock_param(1:irun(5)-1,:);
        ifil(5) = ifil(5) + 1;
        save([DirectoryName '/' M_.fname '_param' int2str(ifil(5)) '.mat'],'stock','stock_logpo','stock_ys');
        if RemoteFlag==1,
            OutputFileName_param = [OutputFileName_param; {[DirectoryName filesep], [M_.fname '_param' int2str(ifil(5)) '.mat']}];
        end
        irun(5) = 1;
    end

    if horizon && (irun(6) > MAX_nforc1 || b == B)
        stock = stock_forcst_mean(:,:,1:irun(6)-1);
        ifil(6) = ifil(6) + 1;
        save([DirectoryName '/' M_.fname '_forc_mean' int2str(ifil(6)) '.mat'],'stock');
        if RemoteFlag==1,
            OutputFileName_forc_mean = [OutputFileName_forc_mean; {[DirectoryName filesep], [M_.fname '_forc_mean' int2str(ifil(6)) '.mat']}];
        end
        irun(6) = 1;
    end

    if horizon && (irun(7) > MAX_nforc2 ||  b == B)
        stock = stock_forcst_point(:,:,1:irun(7)-1);
        ifil(7) = ifil(7) + 1;
        save([DirectoryName '/' M_.fname '_forc_point' int2str(ifil(7)) '.mat'],'stock');
        if RemoteFlag==1,
            OutputFileName_forc_point = [OutputFileName_forc_point; {[DirectoryName filesep], [M_.fname '_forc_point' int2str(ifil(7)) '.mat']}];
        end
        irun(7) = 1;
    end

    % if moments_varendo && (irun(8) > MAX_momentsno || b == B)
    %    stock = stock_moments(1:irun(8)-1);
    %    ifil(8) = ifil(8) + 1;
    %    save([DirectoryName '/' M_.fname '_moments' int2str(ifil(8)) '.mat'],'stock');
    %    if RemoteFlag==1,
    %    OutputFileName_moments = [OutputFileName_moments; {[DirectoryName filesep], [M_.fname '_moments' int2str(ifil(8)) '.mat']}];
    %    end
    %    irun(8) = 1;
    % end

    %   DirectoryName=TempPath;


    if exist('OCTAVE_VERSION'),
        if (whoiam==0),
            printf('Taking subdraws: %3.f%% done\r', b/B*100);
        end
    elseif ~whoiam,
        waitbar(b/B,h);
    end

    if  whoiam,
        if ~exist('OCTAVE_VERSION')
            fprintf('Done! \n');
        end
        waitbarString = [ 'Subdraw ' int2str(b) '/' int2str(B) ' done.'];
        fMessageStatus((b-fpar+1)/(B-fpar+1),whoiam,waitbarString, waitbarTitle, Parallel(ThisMatlab));
    end
end

myoutput.ifil=ifil;
if RemoteFlag==1,
    myoutput.OutputFileName = [OutputFileName_smooth;
                        OutputFileName_update;
                        OutputFileName_inno;
                        OutputFileName_error;
                        OutputFileName_filter_step_ahead;
                        OutputFileName_param;
                        OutputFileName_forc_mean;
                        OutputFileName_forc_point];
    % OutputFileName_moments];
end

if exist('OCTAVE_VERSION')
    printf('\n');
    diary on;
else
    if exist('h')
        close(h)
    end

end
