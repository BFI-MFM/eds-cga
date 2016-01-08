function dynare_estimation(var_list,varargin)
% function dynare_estimation(var_list)
% runs the estimation of the model
%  
% INPUTS
%   var_list:  selected endogenous variables vector
%  
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2003-2011 Dynare Team
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

global options_ oo_ M_ oo_recursive_

%% Decide if a DSGE or DSGE-VAR has to be estimated.
if ~isempty(strmatch('dsge_prior_weight',M_.param_names))
    options_.dsge_var = 1;
end

var_list = check_list_of_variables(options_, M_, var_list);
options_.varlist = var_list;

if isfield(options_,'nobs')
    nobs = options_.nobs;
else
    nobs = [];
end

nnobs = length(nobs);
horizon = options_.forecast;

if nnobs > 1
    if nargin > 1
        dname = vargin{1};
    else
        dname = M_.fname;
    end
    for i=1:nnobs
        options_.nobs = nobs(i);
        dynare_estimation_1(var_list,[dname '_' int2str(nobs(i))]);
        oo_recursive_{nobs(i)} = oo_;
    end
else
    dynare_estimation_1(var_list,varargin{:});
end

if nnobs > 1 && horizon > 0
    mh_replic = options_.mh_replic;
    rawdata = read_variables(options_.datafile,options_.varobs,[],options_.xls_sheet,options_.xls_range);
    gend = options_.nobs;
    rawdata = rawdata(options_.first_obs:options_.first_obs+gend-1,:);
    % Take the log of the variables if needed
    if options_.loglinear && ~options_.logdata   % and if the data are not in logs, then...
        rawdata = log(rawdata);  
    end

    endo_names = M_.endo_names;
    n_varobs = size(options_.varobs,1);
    
    if isempty(var_list)
        var_list = endo_names;
        nvar    = size(endo_names,1);
        SelecVariables = transpose(1:nvar);
    else
        nvar = size(var_list,1);
        SelecVariables = [];
        for i=1:nvar
            if ~isempty(strmatch(var_list(i,:),endo_names,'exact'))
                SelecVariables = [SelecVariables;strmatch(var_list(i,:),endo_names, ...
                                                          'exact')];
            else
                error(['Estimation: ' var_list(i,:) ' isn''t an endogenous' ...
                       'variable'])
            end   
        end
    end

    IdObs    = zeros(n_varobs,1);
    for j=1:n_varobs
        for i=1:nvar
            iobs = strmatch(options_.varobs(j,:),var_list,'exact');
        end
        if ~isempty(iobs)
            IdObs(j,1) = iobs;
        end    
    end 

    k = 3+min(nobs(end)-nobs(1)+horizon, ...
              size(rawdata,1)-nobs(1));
    data2 = rawdata(end-k+1:end,:);
    [nbplt,nr,nc,lr,lc,nstar] = pltorg(nvar);
    m = 1;
    for i = 1:size(var_list,1)
        if mod(i,nstar) == 1
            hfig = figure('Name','Out of sample forecasts');
            m = 1;
        end
        subplot(nr,nc,m)
        hold on
        if any(i==IdObs)
            k2 = find(i==IdObs);
            if options_.loglinear == 1
                plot(1:k,exp(data2(end-k+1:end,k2))','-k','linewidth',2);
            else
                plot(1:k,data2(end-k+1:end,k2)','-k','linewidth',2);
            end
            offsetx = 3;
        else
            offsetx = 0;
        end      
        vname = deblank(var_list(i,:));
        maxlag = M_.maximum_lag;
        for j=1:nnobs
            if mh_replic > 0
                eval(['oo_.RecursiveForecast.Mean.' vname '(j,:) =' ...
                      'oo_recursive_{' int2str(nobs(j))  '}.MeanForecast.Mean.' ...
                      vname '(maxlag+1:end);']);
                eval(['oo_.RecursiveForecast.HPDinf.' vname '(j,:) =' ...
                      'oo_recursive_{' int2str(nobs(j))  '}.MeanForecast.HPDinf.' ...
                      vname '(maxlag+1:end);']);
                eval(['oo_.RecursiveForecast.HPDsup.' vname '(j,:) =' ...
                      'oo_recursive_{' int2str(nobs(j))  '}.MeanForecast.HPDsup.' ...
                      vname '(maxlag+1:end);']);
                eval(['oo_.RecursiveForecast.HPDTotalinf.' vname '(j,:) =' ...
                      'oo_recursive_{' int2str(nobs(j))  '}.PointForecast.HPDinf.' ...
                      vname '(maxlag+1:end);']);
                eval(['oo_.RecursiveForecast.HPDTotalsup.' vname '(j,:) =' ...
                      'oo_recursive_{' int2str(nobs(j))  '}.PointForecast.HPDsup.' ...
                      vname '(maxlag+1:end);']);
            else
                eval(['oo_.RecursiveForecast.Mean.' vname '(j,:) =' ...
                      'oo_recursive_{' int2str(nobs(j))  '}.forecast.Mean.' ...
                      vname ';']);
                eval(['oo_.RecursiveForecast.HPDinf.' vname '(j,:) =' ...
                      'oo_recursive_{' int2str(nobs(j))  '}.forecast.HPDinf.' ...
                      vname ';']);
                eval(['oo_.RecursiveForecast.HPDsup.' vname '(j,:) =' ...
                      'oo_recursive_{' int2str(nobs(j))  '}.forecast.HPDsup.' ...
                      vname ';']);
            end
            x = offsetx+nobs(j)-nobs(1)+(1:horizon);
            y = eval(['oo_.RecursiveForecast.Mean.' vname '(j,:)']);
            y1 = eval(['oo_.RecursiveForecast.HPDinf.' vname '(j,:)']);
            y2 = eval(['oo_.RecursiveForecast.HPDsup.' vname '(j,:)']);
            plot(x,y,'-b','linewidth',2)
            plot(x,y1,'--g', ...
                 'linewidth',1.5)
            plot(x,y2,'--g', ...
                 'linewidth',1.5)
            if mh_replic
                y3 = eval(['oo_.RecursiveForecast.HPDTotalinf.' vname '(j,:)']);
                y4 = eval(['oo_.RecursiveForecast.HPDTotalsup.' vname ...
                           '(j,:)']);
                plot(x,y3,'--r', ...
                     'linewidth',1.5)
                plot(x,y4,'--r','linewidth',1.5)
            end
        end
        %    set(gca,'XTick',offsetx+[1 10 20 30 40 50 60 70 80 90]);
        %    set(gca,'XTickLabel',{'1';'10';'20';'30';'40';'50';'60';'70';'80';'90'});
        %   xlim([1 options_.forecast+10]);
        if any(k==IdObs)
            plot([offsetx+1 offsetx+1],ylim,'-c')
        end
        box on
        title(vname,'Interpreter','none')
        hold off
        m = m + 1;
    end
end
