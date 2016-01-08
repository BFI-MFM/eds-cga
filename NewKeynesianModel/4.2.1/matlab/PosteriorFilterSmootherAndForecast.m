function PosteriorFilterSmootherAndForecast(Y,gend, type,data_index)

% function PosteriorFilterSmootherAndForecast(Y,gend, type)
% Computes posterior filter smoother and forecasts
%
% INPUTS
%    Y:      data
%    gend:   number of observations 
%    type:   posterior
%            prior
%            gsa
%    
% OUTPUTS
%    none
%        
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2005-2010 Dynare Team
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

nvx  = estim_params_.nvx;
nvn  = estim_params_.nvn;
ncx  = estim_params_.ncx;
ncn  = estim_params_.ncn;
np   = estim_params_.np ;
npar = nvx+nvn+ncx+ncn+np;
offset = npar-np;
naK = length(options_.filter_step_ahead);
%%
MaxNumberOfPlotPerFigure = 4;% The square root must be an integer!
MaxNumberOfBytes=options_.MaxNumberOfBytes;
endo_nbr=M_.endo_nbr;
exo_nbr=M_.exo_nbr;
nvobs     = size(options_.varobs,1);
nn = sqrt(MaxNumberOfPlotPerFigure);
iendo = 1:endo_nbr;
i_last_obs = gend+(1-M_.maximum_endo_lag:0);
horizon = options_.forecast;
maxlag = M_.maximum_endo_lag;
%%
CheckPath('Plots/');
DirectoryName = CheckPath('metropolis');
load([ DirectoryName '/'  M_.fname '_mh_history.mat'])
FirstMhFile = record.KeepedDraws.FirstMhFile;
FirstLine = record.KeepedDraws.FirstLine; 
TotalNumberOfMhFiles = sum(record.MhDraws(:,2)); LastMhFile = TotalNumberOfMhFiles; 
TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
clear record;
B = min(1200, round(0.25*NumberOfDraws));
B = 200;
%%
MAX_nruns = min(B,ceil(options_.MaxNumberOfBytes/(npar+2)/8));
MAX_nsmoo = min(B,ceil(MaxNumberOfBytes/((endo_nbr)*gend)/8));
MAX_ninno = min(B,ceil(MaxNumberOfBytes/(exo_nbr*gend)/8));
MAX_nerro = min(B,ceil(MaxNumberOfBytes/(size(options_.varobs,1)*gend)/8));
if naK
    MAX_naK   = min(B,ceil(MaxNumberOfBytes/(size(options_.varobs,1)* ...
                                             length(options_.filter_step_ahead)*gend)/8));
end
if horizon
    MAX_nforc1 = min(B,ceil(MaxNumberOfBytes/((endo_nbr)*(horizon+maxlag))/8));
    MAX_nforc2 = min(B,ceil(MaxNumberOfBytes/((endo_nbr)*(horizon+maxlag))/ ...
                            8));
    IdObs    = bayestopt_.mfys;

end

%%
varlist = options_.varlist;
if isempty(varlist)
    varlist = M_.endo_names;
    SelecVariables = transpose(1:M_.endo_nbr);
    nvar = M_.endo_nbr;
else
    nvar = size(varlist,1);
    SelecVariables = [];
    for i=1:nvar
        if ~isempty(strmatch(varlist(i,:),M_.endo_names,'exact'))
            SelecVariables = [SelecVariables;strmatch(varlist(i,:),M_.endo_names,'exact')];
        end
    end
end

irun1 = 1;
irun2 = 1;
irun3 = 1;
irun4 = 1;
irun5 = 1;
irun6 = 1;
irun7 = 1;
ifil1 = 0;
ifil2 = 0;
ifil3 = 0;
ifil4 = 0;
ifil5 = 0;
ifil6 = 0;
ifil7 = 0;

h = waitbar(0,'Bayesian smoother...');

stock_param = zeros(MAX_nruns, npar);
stock_logpo = zeros(MAX_nruns,1);
stock_ys = zeros(MAX_nruns,endo_nbr);
if options_.smoother
    stock_smooth = zeros(endo_nbr,gend,MAX_nsmoo);
    stock_innov  = zeros(exo_nbr,gend,B);
    stock_error = zeros(nvobs,gend,MAX_nerro);
end
if options_.filter_step_ahead
    stock_filter = zeros(naK,endo_nbr,gend+options_.filter_step_ahead(end),MAX_naK);
end
if options_.forecast
    stock_forcst_mean = zeros(endo_nbr,horizon+maxlag,MAX_nforc1);
    stock_forcst_total = zeros(endo_nbr,horizon+maxlag,MAX_nforc2);
end

for b=1:B
    %deep = GetOneDraw(NumberOfDraws,FirstMhFile,LastMhFile,FirstLine,MAX_nruns,DirectoryName);
    [deep, logpo] = GetOneDraw(type);
    set_all_parameters(deep);
    dr = resol(oo_.steady_state,0);
    [alphahat,etahat,epsilonhat,ahat,SteadyState,trend_coeff,aK] = ...
        DsgeSmoother(deep,gend,Y,data_index);
    
    if options_.loglinear
        stock_smooth(dr.order_var,:,irun1) = alphahat(1:endo_nbr,:)+ ...
            repmat(log(dr.ys(dr.order_var)),1,gend);
    else
        stock_smooth(dr.order_var,:,irun1) = alphahat(1:endo_nbr,:)+ ...
            repmat(dr.ys(dr.order_var),1,gend);
    end    
    if nvx
        stock_innov(:,:,irun2)  = etahat;
    end
    if nvn
        stock_error(:,:,irun3)  = epsilonhat;
    end
    if naK
        stock_filter(:,dr.order_var,:,irun4) = aK(options_.filter_step_ahead,1:endo_nbr,:);
    end
    stock_param(irun5,:) = deep;
    stock_logpo(irun5,1) = logpo;
    stock_ys(irun5,:) = SteadyState';

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
            %      yf = exp(yf);
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
            %      yf1 = exp(yf1);
        else
            yf1 = yf1 + repmat(SteadyState',[horizon+maxlag,1,1]);
        end

        stock_forcst_mean(:,:,irun6) = yf';
        stock_forcst_total(:,:,irun7) = yf1';
    end
    
    irun1 = irun1 + 1;
    irun2 = irun2 + 1;
    irun3 = irun3 + 1;
    irun4 = irun4 + 1;
    irun5 = irun5 + 1;
    irun6 = irun6 + 1;
    irun7 = irun7 + 1;

    if irun1 > MAX_nsmoo | b == B
        stock = stock_smooth(:,:,1:irun1-1);
        ifil1 = ifil1 + 1;
        save([DirectoryName '/' M_.fname '_smooth' int2str(ifil1) '.mat'],'stock');
        irun1 = 1;
    end
    
    if nvx & (irun2 > MAX_ninno | b == B)
        stock = stock_innov(:,:,1:irun2-1);
        ifil2 = ifil2 + 1;
        save([DirectoryName '/' M_.fname '_inno' int2str(ifil2) '.mat'],'stock');
        irun2 = 1;
    end
    
    if nvn & (irun3 > MAX_error | b == B)
        stock = stock_error(:,:,1:irun3-1);
        ifil3 = ifil3 + 1;
        save([DirectoryName '/' M_.fname '_error' int2str(ifil3) '.mat'],'stock');
        irun3 = 1;
    end
    
    if naK & (irun4 > MAX_naK | b == B)
        stock = stock_filter(:,:,:,1:irun4-1);
        ifil4 = ifil4 + 1;
        save([DirectoryName '/' M_.fname '_filter' int2str(ifil4) '.mat'],'stock');
        irun4 = 1;
    end
    
    if irun5 > MAX_nruns | b == B
        stock = stock_param(1:irun5-1,:);
        ifil5 = ifil5 + 1;
        save([DirectoryName '/' M_.fname '_param' int2str(ifil5) '.mat'],'stock','stock_logpo','stock_ys');
        irun5 = 1;
    end

    if horizon & (irun6 > MAX_nforc1 | b == B)
        stock = stock_forcst_mean(:,:,1:irun6-1);
        ifil6 = ifil6 + 1;
        save([DirectoryName '/' M_.fname '_forc_mean' int2str(ifil6) '.mat'],'stock');
        irun6 = 1;
    end

    if horizon & (irun7 > MAX_nforc2 |  b == B)
        stock = stock_forcst_total(:,:,1:irun7-1);
        ifil7 = ifil7 + 1;
        save([DirectoryName '/' M_.fname '_forc_total' int2str(ifil7) '.mat'],'stock');
        irun7 = 1;
    end

    waitbar(b/B,h);
end
close(h)

stock_gend=gend;
stock_data=Y;
save([DirectoryName '/' M_.fname '_data.mat'],'stock_gend','stock_data');

if options_.smoother
    pm3(endo_nbr,gend,ifil1,B,'Smoothed variables',...
        M_.endo_names(SelecVariables),M_.endo_names,'tit_tex',M_.endo_names,...
        'names2','smooth',[M_.fname '/metropolis'],'_smooth')
end

if options_.forecast
    pm3(endo_nbr,horizon+maxlag,ifil6,B,'Forecasted variables (mean)',...
        M_.endo_names(SelecVariables),M_.endo_names,'tit_tex',M_.endo_names,...
        'names2','smooth',[M_.fname '/metropolis'],'_forc_mean')
    pm3(endo_nbr,horizon+maxlag,ifil6,B,'Forecasted variables (total)',...
        M_.endo_names(SelecVariables),M_.endo_names,'tit_tex',M_.endo_names,...
        'names2','smooth',[M_.fname '/metropolis'],'_forc_total')
end
