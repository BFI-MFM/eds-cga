function [result,info] = check
% function result = check
% checks determinacy conditions by computing the eigenvalues
%
% INPUTS
%    none
%    
% OUTPUTS
%    result  [integer]   scalar, equal to 1 if the derterministic steady state satisfies BK conditions. 
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


temp_options = options_;
tempex = oo_.exo_simul;
if ~options_.initval_file && M_.exo_nbr > 1
    oo_.exo_simul = ones(M_.maximum_lead+M_.maximum_lag+1,1)*oo_.exo_steady_state';
end

options_.order = 1;

if isempty(options_.qz_criterium)
    options_.qz_criterium = 1+1e-6;
end

[dr, info] = resol(oo_.steady_state,1);

oo_.dr = dr;

if info(1) ~= 0 && info(1) ~= 3 && info(1) ~= 4
    print_info(info, options_.noprint);
end  

oo_.exo_simul = tempex;

eigenvalues_ = dr.eigval;
if (options_.block)
    nyf = dr.nfwrd+dr.nboth;
else
    nyf = nnz(dr.kstate(:,2)>M_.maximum_endo_lag+1);
end;
[m_lambda,i]=sort(abs(eigenvalues_));
n_explod = nnz(abs(eigenvalues_) > options_.qz_criterium);

result = 0;
if (nyf== n_explod) && (dr.rank == nyf)
    result = 1;
end

if options_.noprint == 0
    disp(' ')
    disp('EIGENVALUES:')
    disp(sprintf('%16s %16s %16s\n','Modulus','Real','Imaginary'))
    z=[m_lambda real(eigenvalues_(i)) imag(eigenvalues_(i))]';
    disp(sprintf('%16.4g %16.4g %16.4g\n',z))
    disp(sprintf('\nThere are %d eigenvalue(s) larger than 1 in modulus ', n_explod));
    disp(sprintf('for %d forward-looking variable(s)',nyf));
    disp(' ')
    if dr.rank == nyf && nyf == n_explod
        disp('The rank condition is verified.')
    else
        disp('The rank conditions ISN''T verified!')
    end
    disp(' ')
end

options_ = temp_options;