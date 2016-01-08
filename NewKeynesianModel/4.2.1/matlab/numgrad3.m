function [g, badg, f0, f1, f2] = numgrad3(fcn,f0,x,epsilon,varargin)
% Computes the gradient of the objective function fcn using a three points
% formula if possible.
%
% Adapted from Sims' numgrad routine.
%
% See section 25.3.4 in Abramovitz and Stegun (1972, Tenth Printing, December) Handbook of Mathematical Functions.
% http://www.math.sfu.ca/~cbm/aands/ 

% Original file downloaded from:
% http://sims.princeton.edu/yftp/optimize/mfiles/numgrad.m

% Copyright (C) 1993-2007 Christopher Sims
% Copyright (C) 2008-2010 Dynare Team
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

f1 = NaN;
f2 = NaN;

delta = epsilon;
n=length(x);
tvec=delta*eye(n);
g=zeros(n,1);

badg=0;
goog=1;
scale=1;
for i=1:n
    if size(x,1)>size(x,2)
        tvecv=tvec(i,:);
    else
        tvecv=tvec(:,i);
    end
    [f1,cost_flag1] = feval(fcn, x+scale*transpose(tvecv), varargin{:});
    [f2,cost_flag2] = feval(fcn, x-scale*transpose(tvecv), varargin{:});
    if cost_flag1 && cost_flag2
        g0 = (f1 - f2) / (2*scale*delta);
    else
        if cost_flag1
            g0 = (f1-f0) / (scale*delta);
        elseif cost_flag2
            g0 = (f0-f2) / (scale*delta);
        else
            goog=0;
        end
    end
    if goog && abs(g0)< 1e15 
        g(i)=g0;
    else
        disp('bad gradient ------------------------')
        g(i)=0;
        badg=1;
    end
end