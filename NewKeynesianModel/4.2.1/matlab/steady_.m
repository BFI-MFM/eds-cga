function steady_()
% function steady_()
% Computes the steady state 
%  
% INPUTS
%   none
%  
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

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

global M_ oo_ options_

if options_.solve_algo < 0 || options_.solve_algo > 8
    error('STEADY: solve_algo must be between 0 and 8')
end

if ~options_.bytecode && ~options_.block && options_.solve_algo > 4
    error('STEADY: you can''t use solve_algo > 4 without block nor bytecode options')
end

if ~options_.bytecode && options_.block && options_.solve_algo == 5
    error('STEADY: you can''t use solve_algo = 5 without bytecode option')
end

if exist('OCTAVE_VERSION') && options_.solve_algo == 7
    error('SIMUL: you can''t use solve_algo = 7 under Octave')
end

if options_.steadystate_flag
    [ys,check] = feval([M_.fname '_steadystate'],...
                       oo_.steady_state,...
                       [oo_.exo_steady_state; ...
                        oo_.exo_det_steady_state]);
    if size(ys,1) < M_.endo_nbr 
        if length(M_.aux_vars) > 0
            ys = add_auxiliary_variables_to_steadystate(ys,M_.aux_vars,...
                                                        M_.fname,...
                                                        oo_.exo_steady_state,...
                                                        oo_.exo_det_steady_state,...
                                                        M_.params,...
                                                        options_.bytecode);
        else
            error([M_.fname '_steadystate.m doesn''t match the model']);
        end
    end
    oo_.steady_state = ys;
    % Check if the steady state obtained from the _steadystate file is a steady state.
    check1 = 0;
    if isempty(options_.unit_root_vars)
        if options_.block && ~options_.bytecode 
            check2 = zeros(size(M_.blocksMFS,1),1);
            for b = 1:size(M_.blocksMFS,1)
                n = size(M_.blocksMFS{b}, 1);
                ss = oo_.steady_state;
                if n
                    check2(b) = max(abs(feval([M_.fname '_static'], b, ss, ...
                                              [oo_.exo_steady_state; oo_.exo_det_steady_state], M_.params))) > options_.dynatol;
                else
                    [r, g1, ssss] = feval([M_.fname '_static'], b, ss, ...
                                          [oo_.exo_steady_state; oo_.exo_det_steady_state], M_.params);
                end
            end
            check1 = any(check2);
            idx = find(abs(ssss-oo_.steady_state)>10*options_.dynatol);
            if ~isempty(idx)
                check1 = 1;
            end
        elseif options_.bytecode
            [check1, residuals] = bytecode('evaluate','static',oo_.steady_state,...
                                           [oo_.exo_steady_state; ...
                                oo_.exo_det_steady_state], M_.params, oo_.steady_state, 1);
            mexErrCheck('bytecode', check1);
            check1 = max(abs(residuals)) > options_.dynatol ;
        else
            check1 = 0;
            check1 = max(abs(feval([M_.fname '_static'],...
                                   oo_.steady_state,...
                                   [oo_.exo_steady_state; ...
                                oo_.exo_det_steady_state], M_.params))) > options_.dynatol ;
        end
    end
    if check1
        resid;
        error(['The steadystate values returned by ' M_.fname ...
               '_steadystate.m don''t solve the static model!' ])
    end
    if ~isempty(options_.steadystate_partial)
        ssvar = options_.steadystate_partial.ssvar;
        nov   = length(ssvar);
        indv  = zeros(nov,1);
        for i = 1:nov
            indv(i) = strmatch(ssvar(i),M_.endo_names,'exact');
        end
        [oo_.steady_state,check] = dynare_solve('restricted_steadystate',...
                                                oo_.steady_state(indv),...
                                                options_.jacobian_flag, ...         
                                                [oo_.exo_steady_state;oo_.exo_det_steady_state],indv);
    end
elseif (options_.bytecode == 0 && options_.block == 0)
    [oo_.steady_state,check] = dynare_solve([M_.fname '_static'],...
                                            oo_.steady_state,...
                                            options_.jacobian_flag, ...     
                                            [oo_.exo_steady_state; ...
                        oo_.exo_det_steady_state], M_.params);
else
    [oo_.steady_state,check] = dynare_solve_block_or_bytecode(oo_.steady_state,...
                                                      [oo_.exo_steady_state; ...
                        oo_.exo_det_steady_state], M_.params);
end

if check ~= 0
    error('STEADY: convergence problems')
end
