function read_data_
% function read_data_
% reads endogenous and exogenous variables from a text file 
% Used by datafile option in simulate
%
% INPUT
%   none
%
% OUTPUT
%   none
%
% SPECIAL REQUIREMENT
%   none

% Copyright (C) 2007-2010 Dynare Team
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

global options_ M_ oo_;
dname= options_.datafile;

if size(oo_.endo_simul,2) < M_.maximum_lag+M_.maximum_lead+options_.periods
    fid = fopen([dname '_endo.dat'],'r');
    names_line = fgetl(fid);
    allVariables = '';
    positions = ones(0);
    while (any(names_line))
        [chopped,names_line] = strtok(names_line);
        if isempty(allVariables)
            allVariables = chopped;
        else
            allVariables = char(allVariables, chopped);
        end
        positions = [positions ; strmatch(chopped,M_.endo_names,'exact')];
    end
    Values=fscanf(fid,'%f',inf);
    Values=reshape(Values,M_.orig_endo_nbr,size(Values,1)/M_.orig_endo_nbr);
    oo_.endo_simul=[Values(positions,:); kron(oo_.steady_state((M_.orig_endo_nbr+1) : M_.endo_nbr , 1) , ones(1 , size(Values, 2)))];
    fclose(fid);
end

if size(oo_.exo_simul,1) < M_.maximum_lag+M_.maximum_lead+options_.periods
    fid = fopen([dname '_exo.dat'],'r');
    names_line = fgetl(fid);
    allVariables = '';
    positions = ones(0);
    while (any(names_line))
        [chopped,names_line] = strtok(names_line);
        if isempty(allVariables)
            allVariables = chopped;
        else
            allVariables = char(allVariables, chopped);
        end
        positions = [positions ; strmatch(chopped,M_.exo_names,'exact')];
    end
    Values=fscanf(fid,'%f',inf);
    Values=reshape(Values,M_.exo_nbr,size(Values,1)/M_.exo_nbr);
    oo_.exo_simul=(Values(positions,:))';
    fclose(fid);
end
%disp([allVariables M_.endo_names]);
%disp(positions);

end