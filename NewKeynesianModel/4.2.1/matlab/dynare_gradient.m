function [F,G] = dynare_gradient(fcn,x,epsilon,varargin)
% Computes the gradient of a function from R^m in R^n.
%
% INPUTS:
%  fcn      [string]  name of the matlab's function.
%  x        [double]  m*1 vector (where the gradient is evaluated).
%  epsilon  [double]  scalar or m*1 vector of steps. 
%
% OUTPUTS: 
%  F        [double]  n*1 vector, evaluation of the function at x.
%  G        [double]  n*m matrix, evaluation of the gradient at x.
%
% OUTPUTS
% 
% Copyright (C) 2010-2011 Dynare Team
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

% Evaluate the function at x.
F = feval(fcn, x, varargin{:});

% (G)Set dimensions.
m = length(x);
n = length(F);

% Initialization of the gradient.
G = NaN(length(F),length(x));

if length(epsilon==1)
    H = epsilon*eye(m);
else
    H = diag(epsilon);
end

% Compute the gradient.
for i=1:m
    if size(x,1)>size(x,2)
        h = H(i,:);
    else
        h = H(:,i);
    end   
    [Fh,flag] = feval(fcn, x+transpose(h), varargin{:});
    if flag
        G(:,i) = (Fh-F)/epsilon;
    else
        [Fh,flag] = feval(fcn, x-transpose(h), varargin{:});
        if flag
            G(:,i) = (F-Fh)/epsilon;
        else
            error('-- Bad gradient --')
        end
    end
end