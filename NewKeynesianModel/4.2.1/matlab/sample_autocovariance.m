function [autocov,autocor] = sample_autocovariance(data,q)
% Computes the autocovariance function associated to a time series.
% 
%
% INPUTS 
%
%   data            [double]       T*1 vector of data.
%   q               [integer]      Order of the autocovariance function. 
%    
% OUTPUTS 
%   autocov         [double]       (q+1)*1 vector, autocovariance function (first scalar is the variance).  
%   autocor         [double]       (q+1)*1 vector, autocorrelation function (first scalar is equal to one).
%
% SPECIAL REQUIREMENTS

% Copyright (C) 2003-2008 Dynare Team
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

autocov = zeros(q+1,1);

demeaned_data = data(:) - mean(data(:));
sample_size = length( data(q+1:end) );
lagged_indices = transpose(0:-1:-q);

for t = 1:sample_size
    tt = t+q;
    autocov = autocov + demeaned_data(tt)*demeaned_data(tt+lagged_indices);
end

autocov = autocov/sample_size ;

if nargout>1
    autocor = autocov / autocov(1);
end