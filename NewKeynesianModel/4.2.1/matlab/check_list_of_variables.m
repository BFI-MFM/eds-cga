function varlist = check_list_of_variables(options_, M_, varlist)
% This function defines, if necessary, the list of endogenous variables
% for which the posterior statistics have to be computed. 
% 
%
% INPUTS 
%
%   options_        [structure]    Dynare structure.
%   M_              [structure]    Dynare structure (related to model definition).
%   varlist         [string]       Array of strings with name of the endogenous variables.
%    
% OUTPUTS 
%   varlist         [string] 
%        
% SPECIAL REQUIREMENTS

% Copyright (C) 2003-2010 Dynare Team
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

msg = 0;
if options_.dsge_var && options_.bayesian_irf
    if ~isempty(varlist)
        for i=1:size(varlist,1)
            idx = strmatch(deblank(varlist(i,:)),options_.varobs,'exact');
            if isempty(idx)
                disp([varlist(i,:) ' is not an observed variable!']);
                msg = 1;
            end
        end
        if size(varlist,1)~=size(options_.varobs)
            msg = 1;
        end
        if msg
            disp(' ')
            disp('Posterior IRFs will be computed for all observed variables.')
            disp(' ')
        end
    end
    varlist = options_.varobs;
    return
end

if isempty(varlist)
    disp(' ')
    disp(['You did not declare endogenous variables after the estimation command.'])
    cas = [];
    if options_.bayesian_irf
        cas = 'Posterior IRFs';
    end
    if options_.moments_varendo
        if isempty(cas)
            cas = 'Posterior moments';
        else
            cas = [ cas , ', posterior moments'];
        end
    end
    if options_.smoother
        if isempty(cas)
            cas = 'Posterior smoothed variables';
        else
            cas = [ cas , ', posterior smoothed variables'];
        end
    end
    if options_.smoother
        if isempty(cas)
            cas = 'Posterior smoothed variables';
        else
            cas = [ cas , ', posterior smoothed variables'];
        end
    end
    if ~isempty(options_.filter_step_ahead)
        if isempty(cas)
            cas = 'Posterior k-step ahead filtered variables';
        else
            cas = [ cas , ', posterior k-step ahead filtered variables'];
        end
    end
    if options_.forecast
        if isempty(cas)
            cas = 'Posterior forecasts';
        else
            cas = [ cas , ' and posterior forecats'];
        end
    end
    if ~isempty(cas)
        string = [ cas , ' will be computed for the ' num2str(M_.endo_nbr)  ' endogenous variables'];
        string = [ string ' of your model, this can be very long....']; 
        format_text(string, 10)
        choice = [];
        while isempty(choice)
            disp(' ')
            disp(' ')
            disp('Choose one of the following options:')
            disp(' ')
            disp(' [1] Consider all the endogenous variables.')
            disp(' [2] Consider all the observed endogenous variables.')
            disp(' [3] Stop Dynare and change the mod file.')
            disp(' ')
            choice = input('options [default is 1] =  ');
            if isempty(choice)
                choice=1;
            end
            if choice==1
                varlist = M_.endo_names(1:M_.orig_endo_nbr, :);
            elseif choice==2
                varlist = options_.varobs;
            elseif choice==3
                varlist = NaN;
            else
                disp('')
                disp('YOU HAVE TO ANSWER 1, 2 or 3!')
                disp('')
            end
        end
    end
    if isnan(varlist)
        edit([M_.fname '.mod'])
    end
    disp('')
end



function format_text(remain, max_number_of_words_per_line)
index = 0;
line_of_text = [];
while ~isempty(remain)
    [token, remain] = strtok(remain);
    index = index+1;
    if isempty(line_of_text)
        line_of_text = token;
    else
        line_of_text = [line_of_text , ' ' , token];
    end
    if index==max_number_of_words_per_line
        disp(line_of_text)
        index = 0;
        line_of_text = [];
    end
end
if index<max_number_of_words_per_line
    disp(line_of_text)
end