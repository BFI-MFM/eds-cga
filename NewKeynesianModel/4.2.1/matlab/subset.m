function jndx = subset();

% Copyright (C) 2006-2009 Dynare Team
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

global options_ estim_params_ M_

ExcludedParamNames = options_.ExcludedParams;
VarObs = options_.varobs;
VarExo = M_.exo_names;
info   = options_.ParamSubSet;

nvx = estim_params_.nvx;
nvn = estim_params_.nvn;
ncx = estim_params_.ncx;
ncn = estim_params_.ncn;
np  = estim_params_.np;
nx  = nvx+nvn+ncx+ncn+np;

if strcmpi(info,'All')
    indx = (1:nx)';
elseif strcmpi(info,'DeepParameters')
    indx = (nvx+nvn+ncx+ncn+1:nx)';
elseif strcmpi(info,'StructuralShocks')
    indx = [(1:nvx),nvx+nvn+1:nvx+nvn+ncx]';
elseif strcmpi(info,'StructuralShocksWithoutCorrelations')
    indx = (1:nvx)';
elseif strcmpi(info,'MeasurementErrors')
    indx = [(nvx+1:nvx+nvn),(nvx+nvn+ncx+1:nvx+nvn+ncx+ncn)]';
elseif strcmpi(info,'MeasurementErrorsWithoutCorrelations')
    indx = (nvx+1:nvx+nvn)';
elseif strcmpi(info,'AllWithoutMeasurementErrors')
    indx = [(1:nvx),nvx+nvn+1:nvx+nvn+ncx,nvx+nvn+ncx+ncn+1:nx]';
elseif strcmpi(info,'None')
    indx = [];
end

if isempty(ExcludedParamNames)
    jndx = indx;
else
    tt = [];
    for i = 1:length(ExcludedParamNames)
        tmp = strmatch(ExcludedParamNames{i},M_.exo_names);
        if ~isempty(tmp) & ( strcmpi(info,'All') | strcmpi(info,'StructuralShocks') | ...
                             strcmpi(info,'StructuralShocksWithoutCorrelations') | ...
                             strcmpi(info,'AllWithoutMeasurementErrors') )
            % The parameter the user wants to exclude is related to the size of the structural innovations.
            if ncx
                disp(['I cannot exclude some of the structural variances if the'])
                disp(['structural innovations are correlated...'])
                error
            end
            tt = [tt;tmp];
        elseif isempty(tmp) & nvn 
            tmp = strmatch(ExcludedParamNames{i},options_.varobs);
            if ~isempty(tmp) & ( strcmpi(info,'All') | strcmpi(info,'MeasurementErrors') | ...
                                 strcmpi(info,'MeasurementErrorsWithoutCorrelations') )
                % The parameter the user wants to exclude is related to the size of the measurement errors variances.
                tmp = nvx+tmp;
                if ncn
                    disp(['I cannot exclude some the measurement error variances if the'])
                    disp(['measurement errors are correlated...'])
                    error
                end
                tt = [tt;tmp];
            end
        else% Excluded parameters are deep parameters...
            tmp = strmatch(ExcludedParamNames{i},M_.param_names(estim_params_.param_vals(:,1),:),'exact');
            if ~isempty(tmp)
                tt = [tt;nvx+nvn+ncx+ncn+tmp];
            else
                disp('The parameter you want to exclude is unknown!')
                error
            end
        end
    end
    jndx = [];
    for i=1:length(indx)
        if ~any(indx(i)==tt)
            jndx = [ jndx ; indx(i) ];
        end
    end
end