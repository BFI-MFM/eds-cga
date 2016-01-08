function [info,description] = check_prior_analysis_data(type,M_)
% Copyright (C) 2009 Dynare Team
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

info = 0;
if nargout>1
    description = '';
end

%% Get informations about prior draws files. 
if ~exist([ M_.dname '/prior/draws'],'dir')
    disp('check_prior_analysis_data:: Can''t find any prior draws file!')
    return
end

prior_draws_info = dir([ M_.dname '/prior/draws/prior_draws*.mat']);
name_of_the_last_prior_draw_file = prior_draws_info(end).name;
date_of_the_last_prior_draw_file = prior_draws_info(end).datenum;

%% Get informations about _posterior_draws files.
if isempty(prior_draws_info)
    info = 1;
    if nargout>1
        description = 'prior_sampler has to be called.';
    end
    return
else
    number_of_last_prior_draws_file = length(prior_draws_info);
    date_of_the_prior_definition = get_date_of_a_file([ M_.dname '/prior/definition.mat']);
    if date_of_the_prior_definition>date_of_the_last_prior_draw_file
        info = 2;
        if nargout>1
            description = 'prior draws files have to be updated.';
        end
        return
    else
        info = 3; % Nothing to do!  
        if nargout>1
            description = 'prior draws files are up to date.';
        end
    end
end

%% Get informations about prior data files.
switch type
  case 'variance'
    generic_prior_data_file_name = 'Prior2ndOrderMoments';
  case 'decomposition'
    generic_prior_data_file_name = 'PriorVarianceDecomposition';
  case 'correlation'
    generic_prior_data_file_name = 'PriorCorrelations';
  case 'conditional decomposition'
    generic_prior_data_file_name = 'PriorConditionalVarianceDecomposition';
  otherwise
    disp(['This feature is not yet implemented!'])
end
CheckPath('prior/moments');
pdfinfo = dir([ M_.dname '/prior/' generic_prior_data_file_name '*']);
if isempty(pdfinfo)
    info = 4;
    if nargout>1
        description = 'prior draws files have to be processed.';
    end
    return
else
    number_of_the_last_prior_data_file = length(pdfinfo);
    name_of_the_last_prior_data_file = pdinfo(end).name;
    pdfdate = pdinfo(end).datenum;
    % /!\ REMARK /!\
    % The user can change the model or the value of a calibrated 
    % parameter without changing the prior. In this case the (prior) 
    % moments should be computed. But this case cannot be detected!!! 
    if pdfdate<date_of_the_last_prior_draw_file
        info = 5; % prior data files have to be updated.
        if nargout>1
            description = 'prior data files have to be updated.';
        end
    else
        info = 6; % Ok (nothing to do ;-)
        if nargout>1
            description = 'prior data files are up to date.';
        end
    end
end