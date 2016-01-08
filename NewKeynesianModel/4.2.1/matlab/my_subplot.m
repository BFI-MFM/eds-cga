function my_subplot(i,imax,irow,icol,fig_title)

% function my_subplot(i,imax,irow,icol,fig_title)
% spreads subplots on several figures according to a maximum number of
% subplots per figure
%
% INPUTS
%   i:          subplot number
%   imax:       total number of subplots
%   irow:       maximum number of rows in a figure
%   icol:       maximum number of columns in a figure
%   fig_title:  title to be repeated on each figure
%
% OUTPUT
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2003-2009 Dynare Team
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

nfig_max = irow*icol;
if imax < nfig_max
    icol = ceil(sqrt(imax));
    irow=icol;
    if (icol-1)*(icol-2) >= imax
        irow = icol-2;
        icol = icol-1;
    elseif (icol)*(icol-2) >= imax
        irow = icol-2;
    elseif icol*(icol-1) >= imax
        irow = icol-1;
    end
end

i1 = mod(i-1,nfig_max);
if i1 == 0
    figure('Name',fig_title);
end

subplot(irow,icol,i1+1);