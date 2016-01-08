function z = resid(junk)
% function z = resid(junk)
%
% Computes static residuals associated with the guess values.
% 
% INPUTS
%    junk:   dummy value for backward compatibility
%    
% OUTPUTS
%    z:      residuals
%        
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2001-2011 Dynare Team
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

global M_ options_ oo_

tagname = 'name';
if nargin && ischar(junk)
    tagname = junk;
end


tags  = M_.equations_tags;
istag = 0;
if length(tags)
    istag = 1;
end



steady_state_old = oo_.steady_state;

% If using a steady state file, initialize oo_.steady_state with that file
if options_.steadystate_flag
    [ys,check] = feval([M_.fname '_steadystate'], ...
                       oo_.steady_state, ...
                       [oo_.exo_steady_state; ...
                        oo_.exo_det_steady_state]);
    if size(ys, 1) < M_.endo_nbr 
        if length(M_.aux_vars) > 0
            ys = add_auxiliary_variables_to_steadystate(ys,M_.aux_vars, ...
                                                        M_.fname, ...
                                                        oo_.exo_steady_state, ...
                                                        oo_.exo_det_steady_state, ...
                                                        M_.params,...
                                                        options_.bytecode);
        else
            error([M_.fname '_steadystate.m doesn''t match the model']);
        end
    end
    oo_.steady_state = ys;
end

% Compute the residuals
if options_.block && ~options_.bytecode
    z = zeros(M_.endo_nbr,1);
    for i = 1:length(M_.blocksMFS)
        [r, g, yy, var_indx] = feval([M_.fname '_static'],...
                                     i,...
                                     oo_.steady_state,...
                                     [oo_.exo_steady_state; ...
                            oo_.exo_det_steady_state], M_.params);
        if isempty(M_.blocksMFS{i})
            idx = var_indx;
        else
            idx = M_.blocksMFS{i};
        end
        z(idx) = r;
    end
elseif options_.bytecode
    [check, z] = bytecode('evaluate','static');
    mexErrCheck('bytecode', check);
else
    z = feval([M_.fname '_static'],...
              oo_.steady_state,...
              [oo_.exo_steady_state; ...
               oo_.exo_det_steady_state], M_.params);
end


% Display the non-zero residuals if no return value
if nargout == 0
    for i = 1:4
        disp(' ')
    end
    ind = [];
    disp('Residuals of the static equations:')
    disp(' ')
    for i=1:M_.orig_endo_nbr
        if abs(z(i)) < options_.dynatol/100
            tmp = 0;
        else
            tmp = z(i);
        end
        if istag
            tg = tags(cell2mat(tags(:,1)) == i,2:3); % all tags for equation i
            ind = strmatch( tagname, cellstr( tg(:,1) ) );
        end
        if ~istag || length(ind) == 0
            disp(['Equation number ' int2str(i) ' : ' num2str(tmp)])
        else
            t1 = tg( ind , 2 );
            s = cell2mat(t1);
            disp( ['Equation number ', int2str(i) ,' : ', num2str(tmp) ,' : ' s])
        end
    end
    for i = 1:2
        disp(' ')
    end
end

oo_.steady_state = steady_state_old;
