function [z,zss]=dyn2vec(s1,s2)
% function [z,zss]=dyn2vec(s1,s2)
% Takes Dynare variables from oo_.endo_simul and copies them into matlab global vectors
%
% INPUTS
%    s1:    subset of variables to be saved
%    s2:    optional parameter, copies Dynare variables (s1) in matlab ones (s2)
%
% OUTPUTS
%    z:     subset of oo_.endo_simul
%    zss:   matrix of variables in steady state
%
% SPECIAL REQUIREMENTS
%   none
%  

% Copyright (C) 2001-2009 Dynare Team
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

global M_ oo_ options_

if options_.smpl == 0
    k = [1:size(oo_.endo_simul,2)];
else
    k = [M_.maximum_lag+options_.smpl(1):M_.maximum_lag+options_.smpl(2)];
end

if nargin == 0
    if nargout > 0
        t = ['DYNARE dyn2vec error: the function doesn''t return values when' ...
             ' used without input argument'];
        error(t);
    end
    for i=1:size(oo_.endo_simul,1)
        assignin('base',deblank(M_.endo_names(i,:)),oo_.endo_simul(i,k)');
    end
    return
else
    j = strmatch(s1,M_.endo_names,'exact'); 
    if ~ isempty(j)
        z = oo_.endo_simul(j,k)';
    else
        j = strmatch(s1,M_.exo_names,'exact');
        if ~ isempty(j)
            if options_.smpl == 0
                z = oo_.exo_simul(:,j);
            else
                z = oo_.exo_simul(M_.maximum_lag+options_.smpl(1):M_.maximum_lag+options_.smpl(2));
            end
        else
            t = ['DYNARE dyn2vec error: variable ' deblank(s1(i,:)) ' doesn''t' ...
                 ' exist.'] ;
            error (t) ;
        end
    end
end

if nargout == 0
    if nargin == 1
        assignin('base',s1,z);
    elseif nargin == 2
        assignin('base',s2,z);
    end
else
    zss=oo_.steady_state(j);
end

% 02/23/01 MJ redone, incorporating FC's improvements
% 08/24/01 MJ replaced globlize by internal assignin
% 08/24/01 MJ added 'exact' to strmatch (thanks to David Vavra)
% 01/31/03 MJ added provision for alternative name of variable



