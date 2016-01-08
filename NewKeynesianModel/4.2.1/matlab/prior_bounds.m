function bounds = prior_bounds(bayestopt)
% function bounds = prior_bounds(bayestopt)
% computes bounds for prior density.
%
% INPUTS
%    bayestopt  [structure]  characterizing priors (shape, mean, p1..p4)
%    
% OUTPUTS
%    bounds     [double]      matrix specifying prior bounds (row= parameter, column=upper&lower bound)
%    
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2011 Dynare Team
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

global options_

pshape = bayestopt.pshape;
p3 = bayestopt.p3;
p4 = bayestopt.p4;
p6 = bayestopt.p6;
p7 = bayestopt.p7;
prior_trunc = options_.prior_trunc;

bounds = zeros(length(p6),2);

for i=1:length(p6)
    switch pshape(i)
      case 1
        if prior_trunc == 0
            bounds(i,1) = p3(i);
            bounds(i,2) = p4(i);
        else
            bounds(i,1) = betainv(prior_trunc,p6(i),p7(i))*(p4(i)-p3(i))+p3(i);
            bounds(i,2) = betainv(1-prior_trunc,p6(i),p7(i))* ...
                (p4(i)-p3(i))+p3(i);
        end
      case 2
        if prior_trunc == 0
            bounds(i,1) = p3(i);
            bounds(i,2) = Inf;
        else
            try
                bounds(i,1) = gaminv(prior_trunc,p6(i),p7(i))+p3(i);
                bounds(i,2) = gaminv(1-prior_trunc,p6(i),p7(i))+p3(i);
            catch
                % Workaround for ticket #161
                if exist('OCTAVE_VERSION')
                    error(['Due to a bug in Octave, you must choose other values for mean and/or variance of your prior on ' bayestopt.name{i} ', or use another shape'])
                else
                    rethrow(lasterror)
                end
            end
        end
      case 3
        if prior_trunc == 0
            bounds(i,1) = -Inf;
            bounds(i,2) = Inf;
        else
            bounds(i,1) = norminv(prior_trunc,p6(i),p7(i));
            bounds(i,2) = norminv(1-prior_trunc,p6(i),p7(i));
        end
      case 4
        if prior_trunc == 0
            bounds(i,1) = p3(i);
            bounds(i,2) = Inf;
        else
            try
                bounds(i,1) = 1/sqrt(gaminv(1-prior_trunc, p7(i)/2, 2/p6(i)))+p3(i);
                bounds(i,2) = 1/sqrt(gaminv(prior_trunc, p7(i)/2, ...
                                            2/p6(i)))+p3(i);
            catch
                % Workaround for ticket #161
                if exist('OCTAVE_VERSION')
                    error(['Due to a bug in Octave, you must choose other values for mean and/or variance of your prior on ' bayestopt.name{i} ', or use another shape'])
                else
                    rethrow(lasterror)
                end
            end
        end
      case 5
        if prior_trunc == 0
            bounds(i,1) = p6(i);
            bounds(i,2) = p7(i);
        else
            bounds(i,1) = p6(i)+(p7(i)-p6(i))*prior_trunc;
            bounds(i,2) = p7(i)-(p7(i)-p6(i))*prior_trunc;
        end
      case 6
        if prior_trunc == 0
            bounds(i,1) = p3(i);
            bounds(i,2) = Inf;
        else
            try
                bounds(i,1) = 1/gaminv(1-prior_trunc, p7(i)/2, 2/p6(i))+p3(i);
                bounds(i,2) = 1/gaminv(prior_trunc, p7(i)/2, 2/p6(i))+ p3(i);
            catch
                % Workaround for ticket #161
                if exist('OCTAVE_VERSION')
                    error(['Due to a bug in Octave, you must choose other values for mean and/or variance of your prior on ' bayestopt.name{i} ', or use another shape'])
                else
                    rethrow(lasterror)
                end
            end
        end
      otherwise
        error(sprintf('prior_bounds: unknown distribution shape (index %d, type %d)', i, pshape(i)));
    end
end