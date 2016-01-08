function [g, badg, f0, f1, f2, f3, f4] = numgrad5(fcn,f0,x,epsilon,varargin)
% Computes the gradient of the objective function fcn using a five points
% formula if possible.
%
% Adapted from Sims' numgrad.m routine.
%
% See section 25.3.6 Abramovitz and Stegun (1972, Tenth Printing, December) Handbook of Mathematical Functions.
% http://www.math.sfu.ca/~cbm/aands/ 
%
% TODO Try Four points formula when cost_flag3=0 or cost_flag4=0.

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
f3 = NaN;
f4 = NaN;

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
    if cost_flag1==0 || cost_flag2==0
        cost_flag3 = 0;
        cost_flag4 = 0;
        disp('numgrad:: I cannot use the five points formula!!')
    else
        [f3,cost_flag3] = feval(fcn, x+2*scale*transpose(tvecv), varargin{:});
        [f4,cost_flag4] = feval(fcn, x-2*scale*transpose(tvecv), varargin{:});
    end
    if cost_flag1 && cost_flag2 && cost_flag3 && cost_flag4% Five Points formula
        g0 = (8*(f1 - f2)+ f4-f3) / (12*scale*delta);
    elseif cost_flag3==0 || cost_flag4==0
        if cost_flag1 && cost_flag2% Three points formula
            g0 = (f1-f2)/(2*scale*delta);
        else
            if cost_flag1% Two points formula
                g0 = (f1-f0) / (scale*delta);
            elseif cost_flag2% Two points formula
                g0 = (f0-f2) / (scale*delta);
            else% Bad gradient!
                goog=0;
            end
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