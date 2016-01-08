function myoutput=PosteriorIRF_core1(myinputs,fpar,npar,whoiam, ThisMatlab)
%   PARALLEL CONTEXT
%   This function perfom in parallel a portion of  PosteriorIRF.m code.
%   This is a special kind of parallel function. Unlike of other parallel functions,
%   that running in parallel a 'for' cycle, this function run in parallel a
%   'while' loop! The parallelization of 'while' loop (when possible) is a more
%   sophisticated procedure.
%
%   See also the comment in random_walk_metropolis_hastings_core.m funtion.
%
% INPUTS
%   See the comment in random_walk_metropolis_hastings_core.m funtion.
%
% OUTPUTS
% o myoutput  [struc]
%  Contained:
%  OutputFileName_dsge, OutputFileName_param and OutputFileName_bvardsge.
%
% ALGORITHM
%   Portion of PosteriorIRF.m function. Specifically the 'while' cycle.
%
% SPECIAL REQUIREMENTS.
%   None.
%
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


global options_ estim_params_ oo_ M_ bayestopt_

if nargin<4,
    whoiam=0;
end

% Reshape 'myinputs' for local computation.
% In order to avoid confusion in the name space, the instruction struct2local(myinputs) is replaced by:

IRUN = myinputs.IRUN;
irun =myinputs.irun;
irun2=myinputs.irun2;
nosaddle=myinputs.nosaddle;
type=myinputs.type;
if ~strcmpi(type,'prior'),
    x=myinputs.x;
end

if options_.dsge_var
    gend=myinputs.gend;
    nvobs=myinputs.nvobs;
    NumberOfParametersPerEquation = myinputs.NumberOfParametersPerEquation;
    NumberOfLags = myinputs.NumberOfLags;
    NumberOfLagsTimesNvobs = myinputs.NumberOfLagsTimesNvobs;
    Companion_matrix = myinputs.Companion_matrix;
end

nvar=myinputs.nvar;
IndxVariables=myinputs.IndxVariables;
MAX_nirfs_dsgevar=myinputs.MAX_nirfs_dsgevar;
MAX_nirfs_dsge=myinputs.MAX_nirfs_dsge;
MAX_nruns=myinputs.MAX_nruns;

NumberOfIRFfiles_dsge=myinputs.NumberOfIRFfiles_dsge;
NumberOfIRFfiles_dsgevar=myinputs.NumberOfIRFfiles_dsgevar;
ifil2=myinputs.ifil2;

if whoiam
    Parallel=myinputs.Parallel;
end


% MhDirectoryName = myinputs.MhDirectoryName;
if strcmpi(type,'posterior')
    MhDirectoryName = CheckPath('metropolis');
elseif strcmpi(type,'gsa')
    if options_.opt_gsa.pprior
        MhDirectoryName = CheckPath(['GSA' filesep 'prior']);
    else
        MhDirectoryName = CheckPath(['GSA' filesep 'mc']);
    end
else
    MhDirectoryName = CheckPath('prior');
end

RemoteFlag = 0;

if whoiam
    waitbarString = ['Please wait... Bayesian (posterior) IRFs computing. (' int2str(fpar) 'of' int2str(npar) ')...'];
    if Parallel(ThisMatlab).Local,
        waitbarTitle=['Local '];
    else
        waitbarTitle=[Parallel(ThisMatlab).ComputerName];
        RemoteFlag =1;
    end
    fMessageStatus(0,whoiam,waitbarString, waitbarTitle, Parallel(ThisMatlab));
else
    if exist('OCTAVE_VERSION')
        diary off;
        printf('\n')
    else
        if strcmpi(type,'posterior')
            h = waitbar(0,'Bayesian (posterior) IRFs...');
            set(h,'Name','Bayesian (posterior) IRFs.');
        elseif strcmpi(type,'gsa')
            h = waitbar(0,'GSA (prior) IRFs...');
        else
            h = waitbar(0,'Bayesian (prior) IRFs...');
        end
    end

end

OutputFileName_bvardsge = {};
OutputFileName_dsge = {};
OutputFileName_param = {};


fpar0=fpar;
fpar = fpar-1;

if whoiam
    ifil2=ifil2(whoiam);
    NumberOfIRFfiles_dsge=NumberOfIRFfiles_dsge(whoiam);
    NumberOfIRFfiles_dsgevar=NumberOfIRFfiles_dsgevar(whoiam);
end

% Parallel 'while' very good!!!

while fpar<npar
    fpar = fpar + 1;
    irun = irun+1;
    irun2 = irun2+1;
    if strcmpi(type,'prior')
        deep = GetOneDraw(type);
    else
        deep = x(fpar,:);
    end
    stock_param(irun2,:) = deep;
    set_parameters(deep);
    [dr,info] = resol(oo_.steady_state,0);
    if info(1)
        nosaddle = nosaddle + 1;
        fpar = fpar - 1;
        irun = irun-1;
        irun2 = irun2-1;
        if info(1) == 1
            errordef = 'Static variables are not uniquely defined';
        elseif info(1) == 2
            errordef = 'Dll problem';
        elseif info(1) == 3
            errordef = 'No stable trajectory';
        elseif info(1) == 4
            errordef = 'Indeterminacy';
        elseif info(1) == 5
            errordef = 'Rank condition  is not satisfied';
        end
        disp(['PosteriorIRF :: Dynare is unable to solve the model (' errordef ')'])
        continue
    end
    SS(M_.exo_names_orig_ord,M_.exo_names_orig_ord) = M_.Sigma_e+1e-14*eye(M_.exo_nbr);
    SS = transpose(chol(SS));
    for i = 1:M_.exo_nbr
        if SS(i,i) > 1e-13
            y=irf(dr,SS(M_.exo_names_orig_ord,i), options_.irf, options_.drop,options_.replic,options_.order);
            if options_.relative_irf
                y = 100*y/cs(i,i);
            end
            for j = 1:nvar
                if max(y(IndxVariables(j),:)) - min(y(IndxVariables(j),:)) > 1e-12
                    stock_irf_dsge(:,j,i,irun) = transpose(y(IndxVariables(j),:));
                end
            end
        end
    end
    if MAX_nirfs_dsgevar
        IRUN = IRUN+1;
        [fval,cost_flag,info,PHI,SIGMAu,iXX] =  DsgeVarLikelihood(deep',gend);
        dsge_prior_weight = M_.params(strmatch('dsge_prior_weight',M_.param_names));
        DSGE_PRIOR_WEIGHT = floor(gend*(1+dsge_prior_weight));
        SIGMA_inv_upper_chol = chol(inv(SIGMAu*gend*(dsge_prior_weight+1)));
        explosive_var  = 1;
        while explosive_var
            % draw from the marginal posterior of SIGMA
            SIGMAu_draw = rand_inverse_wishart(nvobs, DSGE_PRIOR_WEIGHT-NumberOfParametersPerEquation, ...
                                               SIGMA_inv_upper_chol);
            % draw from the conditional posterior of PHI
            PHI_draw = rand_matrix_normal(NumberOfParametersPerEquation,nvobs, PHI, ...
                                          chol(SIGMAu_draw)', chol(iXX)');
            Companion_matrix(1:nvobs,:) = transpose(PHI_draw(1:NumberOfLagsTimesNvobs,:));
            % Check for stationarity
            explosive_var = any(abs(eig(Companion_matrix))>1.000000001);
        end
        % Get the mean
        mu = zeros(1,nvobs);
        % Get rotation
        if dsge_prior_weight > 0
            Atheta(oo_.dr.order_var,M_.exo_names_orig_ord) = oo_.dr.ghu*sqrt(M_.Sigma_e);
            A0 = Atheta(bayestopt_.mfys,:);
            [OMEGAstar,SIGMAtr] = qr2(A0');
        end
        SIGMAu_chol = chol(SIGMAu_draw)';
        SIGMAtrOMEGA = SIGMAu_chol*OMEGAstar';
        PHIpower = eye(NumberOfLagsTimesNvobs);
        irfs = zeros (options_.irf,nvobs*M_.exo_nbr);
        tmp3 = PHIpower(1:nvobs,1:nvobs)*SIGMAtrOMEGA;
        irfs(1,:) = tmp3(:)';
        for t = 2:options_.irf
            PHIpower = Companion_matrix*PHIpower;
            tmp3 = PHIpower(1:nvobs,1:nvobs)*SIGMAtrOMEGA;
            irfs(t,:)  = tmp3(:)'+kron(ones(1,M_.exo_nbr),mu);
        end
        tmp_dsgevar = kron(ones(options_.irf,1),mu);
        for j = 1:(nvobs*M_.exo_nbr)
            if max(irfs(:,j)) - min(irfs(:,j)) > 1e-10
                tmp_dsgevar(:,j) = (irfs(:,j));
            end
        end
        if IRUN < MAX_nirfs_dsgevar
            stock_irf_bvardsge(:,:,:,IRUN) = reshape(tmp_dsgevar,options_.irf,nvobs,M_.exo_nbr);
        else
            stock_irf_bvardsge(:,:,:,IRUN) = reshape(tmp_dsgevar,options_.irf,nvobs,M_.exo_nbr);
            instr = [MhDirectoryName '/' M_.fname '_irf_bvardsge' ...
                     int2str(NumberOfIRFfiles_dsgevar) '.mat stock_irf_bvardsge;'];,
            eval(['save ' instr]);
            if RemoteFlag==1,
                OutputFileName_bvardsge = [OutputFileName_bvardsge; {[MhDirectoryName filesep], [M_.fname '_irf_bvardsge' int2str(NumberOfIRFfiles_dsgevar) '.mat']}];
            end
            NumberOfIRFfiles_dsgevar = NumberOfIRFfiles_dsgevar+1;
            IRUN =0;
            stock_irf_dsgevar = zeros(options_.irf,nvobs,M_.exo_nbr,MAX_nirfs_dsgevar);
        end
    end
    if irun == MAX_nirfs_dsge | irun == npar | fpar == npar
        if fpar == npar
            stock_irf_dsge = stock_irf_dsge(:,:,:,1:irun);
            if MAX_nirfs_dsgevar & (fpar == npar | IRUN == npar)
                stock_irf_bvardsge = stock_irf_bvardsge(:,:,:,1:IRUN);
                instr = [MhDirectoryName '/' M_.fname '_irf_bvardsge' ...
                         int2str(NumberOfIRFfiles_dsgevar) '.mat stock_irf_bvardsge;'];,
                eval(['save ' instr]);
                NumberOfIRFfiles_dsgevar = NumberOfIRFfiles_dsgevar+1;
                if RemoteFlag==1,
                    OutputFileName_bvardsge = [OutputFileName_bvardsge; {[MhDirectoryName filesep], [M_.fname '_irf_bvardsge' int2str(NumberOfIRFfiles_dsgevar) '.mat']}];
                end
                irun = 0;
            end
        end
        save([MhDirectoryName '/' M_.fname '_irf_dsge' int2str(NumberOfIRFfiles_dsge) '.mat'],'stock_irf_dsge');
        if RemoteFlag==1,
            OutputFileName_dsge = [OutputFileName_dsge; {[MhDirectoryName filesep], [M_.fname '_irf_dsge' int2str(NumberOfIRFfiles_dsge) '.mat']}];
        end
        NumberOfIRFfiles_dsge = NumberOfIRFfiles_dsge+1;
        irun = 0;
    end
    if irun2 == MAX_nruns | fpar == npar
        if fpar == npar
            stock_param = stock_param(1:irun2,:);
        end
        stock = stock_param;
        save([MhDirectoryName '/' M_.fname '_param_irf' int2str(ifil2) '.mat'],'stock');
        if RemoteFlag==1,
            OutputFileName_param = [OutputFileName_param; {[MhDirectoryName filesep], [M_.fname '_param_irf' int2str(ifil2) '.mat']}];
        end
        ifil2 = ifil2 + 1;
        irun2 = 0;
    end
    if exist('OCTAVE_VERSION'),
        if (whoiam==0),
            printf(['Posterior IRF  %3.f%% done\r'],(fpar/npar*100));
        end
    elseif ~whoiam,
        waitbar(fpar/npar,h);
    end
    if whoiam,
        if ~exist('OCTAVE_VERSION')
            fprintf('Done! \n');
        end
        waitbarString = [ 'Subdraw ' int2str(fpar) '/' int2str(npar) ' done.'];
        fMessageStatus((fpar-fpar0)/(npar-fpar0),whoiam,waitbarString, waitbarTitle, Parallel(ThisMatlab));
    end
end

if whoiam==0
    if nosaddle
        disp(['PosteriorIRF :: Percentage of discarded posterior draws = ' num2str(nosaddle/(npar+nosaddle))])
    end
    if exist('h')
        close(h);
    end
    if exist('OCTAVE_VERSION')
        printf('\n');
        diary on;
    end
end

% Copy the rusults of computation on the call machine (specifically in the
% directory on call machine that contain the model).

myoutput.OutputFileName = [OutputFileName_dsge;
                    OutputFileName_param;
                    OutputFileName_bvardsge];




