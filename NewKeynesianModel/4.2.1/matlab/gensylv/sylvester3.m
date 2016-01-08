function x=sylvester3(a,b,c,d)
% solves a*x+b*x*c=d

% Copyright (C) 2005-2009 Dynare Team
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

n = size(a,1);
m = size(c,1);
if n == 1
    x=d./(a*ones(1,m)+b*c);
    return
end
if m == 1
    x = (a+c*b)\d;
    return;
end
x=zeros(n,m);
[u,t]=schur(c);
if exist('OCTAVE_VERSION')
    [aa,bb,qq,zz]=qz(full(a),full(b));
    d=qq'*d*u;
else
    [aa,bb,qq,zz]=qz(full(a),full(b),'real'); % available in Matlab version 6.0
    d=qq*d*u;
end
i = 1;
while i < m
    if t(i+1,i) == 0
        if i == 1
            c = zeros(n,1);
        else
            c = bb*(x(:,1:i-1)*t(1:i-1,i));
        end
        x(:,i)=(aa+bb*t(i,i))\(d(:,i)-c);
        i = i+1;
    else
        if i == n
            c = zeros(n,1);
            c1 = zeros(n,1);
        else
            c = bb*(x(:,1:i-1)*t(1:i-1,i));
            c1 = bb*(x(:,1:i-1)*t(1:i-1,i+1));
        end
        z = [aa+bb*t(i,i) bb*t(i+1,i); bb*t(i,i+1) aa+bb*t(i+1,i+1)]...
            \[d(:,i)-c;d(:,i+1)-c1];
        x(:,i) = z(1:n);
        x(:,i+1) = z(n+1:end);
        i = i + 2;
    end
end
if i == m
    c = bb*(x(:,1:m-1)*t(1:m-1,m));
    x(:,m)=(aa+bb*t(m,m))\(d(:,m)-c);
end
x=zz*x*u';

% 01/25/03 MJ corrected bug for i==m (sign of c in x determination)
% 01/31/03 MJ added 'real' to qz call
