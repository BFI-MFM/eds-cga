function print_info(info,noprint)
% Prints error messages
%
% INPUTS
%   info    [double]   vector returned by resol.m 
%   noprint [integer]  equal to 0 if the error message has to be printed. 
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2005-2011 Dynare Team
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

if ~noprint
    switch info(1)
      case 1
        error(['The model doesn''t determine the current variables' ...
               ' uniquely'])
      case 2
        error(['MJDGGES returns the following error code: ' ...
               int2str(info(2))])
      case 3
        error(['Blanchard Kahn conditions are not satisfied: no stable' ...
               ' equilibrium'])
      case 4
        error(['Blanchard Kahn conditions are not satisfied:' ...
               ' indeterminacy'])
      case 5
        error(['Blanchard Kahn conditions are not satisfied:' ...
               ' indeterminacy due to rank failure'])
      case 6
        error('The jacobian matrix evaluated at the steady state is complex')
      case 19
        error('The steadystate file did not compute the steady state (inconsistent deep parameters).')
      case 20
        error(['Impossible to find the steady state. Either the model' ...
               ' doesn''t have a unique steady state of the guess values' ...
               ' are too far from the solution'])
      case 21
        error('The steady state is complex.')
      case 30 
        error('Variance can''t be computed')
      case 41
        error('one (many) parameter(s) do(es) not satisfy the lower bound');
      case 42
        error('one (many) parameter(s) do(es) not satisfy the upper bound');
      case 43
        error('Covariance matrix of shocks is not positive definite')
      case 44 %DsgeLikelihood_hh / DsgeLikelihood
        error('');
      case 51
        error('You are estimating a DSGE-VAR model, but the value of the dsge prior weight is too low!')
      case 52 %DsgeVarLikelihood
        error('');

        % Aim Code Conversions by convertAimCodeToInfo.m
      case 102
        error('Aim: roots not correctly computed by real_schur.');
      case 103
        error('Aim: too many big roots.');
      case 135
        error('Aim: too many big roots, and q(:,right) is singular.');
      case 104
        error('Aim: too few big roots.');
      case 145
        error('Aim: too few big roots, and q(:,right) is singular.');
      case 105
        error('Aim: q(:,right) is singular.');
      case 161
        error('Aim: too many exact shiftrights.');
      case 162
        error('Aim: too many numeric shiftrights.');
      otherwise
        error('This case shouldn''t happen. Contact the authors of Dynare')
    end
end