function [d1,d2] = get_date_of_a_file(filename)
%function [d1,d2] = get_date_of_a_file(filename)

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

info = dir(filename);
if isempty(info)
    error(['get_date_of_a_file:: I''m not able to find ' filename '!'])
end
d1 = info.datenum;
if nargout>1
    d2 = info.date;
end