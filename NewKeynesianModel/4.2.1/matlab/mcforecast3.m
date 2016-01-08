function forcs = mcforecast3(cL,H,mcValue,shocks,forcs,T,R,mv,mu)

% Copyright (C) 2006-2008 Dynare Team
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

if cL
    e = zeros(size(mcValue,1),cL);
    for t = 1:cL
        e(:,t) = inv(mv*R*mu)*(mcValue(:,t)-mv*T*forcs(:,t)-mv*R*shocks(:,t));
        forcs(:,t+1) = T*forcs(:,t)+R*(mu*e(:,t)+shocks(:,t));
    end
end
for t = cL+1:H
    forcs(:,t+1) = T*forcs(:,t)+R*shocks(:,t);
end