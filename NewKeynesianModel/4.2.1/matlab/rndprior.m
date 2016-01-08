function y = rndprior(bayestopt_)
% function y = rndprior(bayestopt_)
% Draws random number from the prior density
%
% INPUTS
%   bayestopt_:    structure characterizing priors
%    
% OUTPUTS
%   y:             drawn numbers vector              
%
% SPECIAL REQUIREMENTS
%    none

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

pshape=bayestopt_.pshape;
p3=bayestopt_.p3;
p4=bayestopt_.p4;
p6=bayestopt_.p6;
p7=bayestopt_.p7;

y = NaN(1,length(pshape)); 

for i=1:length(pshape)
    switch pshape(i)
      case 1 % Beta
        y(i) = betarnd(p6(i),p7(i));
        y(i) = y(1,i) * (p4(i)-p3(i)) + p3(i);
      case 2 % Generalized gamma
        y(i) = gamrnd(p6(i),p7(i)) + p3(i);
      case 3 % Gaussian
        y(i) = randn*p7(i) + p6(i) ;
      case 4 % Inverse-gamma type 1
        y(i) = 1/sqrt(gamrnd(p7(i)/2, 2/p6(i))) + p3(i);
      case 5 % Uniform
        y(i) = rand*(p4(i)-p3(i)) + p3(i);
      case 6 % Inverse-gamma type 2
        y(i) = 1/gamrnd(p7(i)/2, 2/p6(i)) + p3(i); 
      otherwise
        error(sprintf('rndprior: unknown distribution shape (index %d, type %d)', i, pshape(i)));
    end
end