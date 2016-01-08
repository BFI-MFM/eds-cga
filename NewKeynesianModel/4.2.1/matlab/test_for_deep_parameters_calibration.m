function test_for_deep_parameters_calibration(M_)
% Issues a warning is some of the parameters are NaNs.
%
% INPUTS
%   M_    [structure]   Description of the (simulated or estimated) model.
%  
% OUTPUTS
%   none
%    
% ALGORITHM
%   none
%    
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2010 Dynare Team
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
plist = list_of_parameters_calibrated_as_NaN(M_);
if ~isempty(plist)
    message = ['Some of the parameters have no value (' ];
    for i=1:size(plist,1)
        if i<size(plist,1)
            message = [message, deblank(plist(i,:)) ', '];
        else
            message = [message, deblank(plist(i,:)) ')'];
        end
    end
    tmp = dbstack;
    message = [message, ' when using ' tmp(2).name '. '];
    message = [message, 'If these parameters are not initialized in a steadystate file, Dynare may not be able to solve the model...'];
    message_id  = 'Dynare:ParameterCalibration:NaNValues';
    warning(message_id,message);
end