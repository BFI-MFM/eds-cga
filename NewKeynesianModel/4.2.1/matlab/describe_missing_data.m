function [data_index,number_of_observations,no_more_missing_observations] = describe_missing_data(data,gend,nvarobs)

% Copyright (C) 2008-2009 Dynare Team
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

[variable_index,observation_index] = find(~isnan(data));

data_index = cell(1,gend);
missing_observations_counter = NaN(gend,1);
for obs=1:gend
    idx = find(observation_index==obs);
    tmp = variable_index(idx);
    missing_observations_counter(obs,1) = nvarobs-length(tmp);
    data_index(obs) = { tmp(:) };
end
missing_observations_counter = cumsum(missing_observations_counter);

number_of_observations = length(variable_index);

if ~missing_observations_counter
    no_more_missing_observations = 0;
else
    tmp = find(missing_observations_counter>=(gend*nvarobs-number_of_observations));
    no_more_missing_observations = tmp(1);
end