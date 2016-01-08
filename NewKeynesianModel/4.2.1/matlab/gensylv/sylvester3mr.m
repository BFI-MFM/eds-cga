function x=sylvester3mr(a,b,c,d)
% solves a*x+b*x*c=d where d is [n x m x p]

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
if length(size(d))==2,
    x=sylvester3(a,b,c,d);
    return
end
p = size(d,3);
if n == 1
    for j=1:p,
        x(:,:,j)=d(:,:,j)./(a*ones(1,m)+b*c);
    end
    return
end
if m == 1
    invacb = inv(a+c*b);
    x = invacb*d;
    return;
end
x=zeros(n,m,p);
[u,t]=schur(c);
if exist('OCTAVE_VERSION')
    [aa,bb,qq,zz]=qz(full(a),full(b));
    for j=1:p,
        d(:,:,j)=qq'*d(:,:,j)*u;
    end
else
    [aa,bb,qq,zz]=qz(full(a),full(b),'real'); % available in Matlab version 6.0
    for j=1:p,
        d(:,:,j)=qq*d(:,:,j)*u;
    end
end
i = 1;
c = zeros(n,1,p);
while i < m
    if t(i+1,i) == 0
        if i == 1
            c = zeros(n,1,p);
        else
            for j=1:p,
                c(:,:,j) = bb*(x(:,1:i-1,j)*t(1:i-1,i));
            end
        end
        %       aabbtinv = inv(aa+bb*t(i,i));
        %       x(:,i,:)=aabbtinv*squeeze(d(:,i,:)-c);
        x(:,i,:)=(aa+bb*t(i,i))\squeeze(d(:,i,:)-c);
        i = i+1;
    else
        if i == n
            c = zeros(n,1,p);
            c1 = zeros(n,1,p);
        else
            for j=1:p,
                c(:,:,j) = bb*(x(:,1:i-1,j)*t(1:i-1,i));
                c1(:,:,j) = bb*(x(:,1:i-1,j)*t(1:i-1,i+1));
            end
        end
        %       bigmatinv = inv([aa+bb*t(i,i) bb*t(i+1,i); bb*t(i,i+1) aa+bb*t(i+1,i+1)]);
        %       z = bigmatinv * squeeze([d(:,i,:)-c;d(:,i+1,:)-c1]);
        bigmat = ([aa+bb*t(i,i) bb*t(i+1,i); bb*t(i,i+1) aa+bb*t(i+1,i+1)]);
        z = bigmat\squeeze([d(:,i,:)-c;d(:,i+1,:)-c1]);
        x(:,i,:) = z(1:n,:);
        x(:,i+1,:) = z(n+1:end,:);
        i = i + 2;
    end
end
if i == m
    for j=1:p,
        c(:,:,j) = bb*(x(:,1:m-1,j)*t(1:m-1,m));
    end
    %       aabbtinv = inv(aa+bb*t(m,m));
    %     x(:,m,:)=aabbtinv*squeeze(d(:,m,:)-c);
    aabbt = (aa+bb*t(m,m));
    x(:,m,:)=aabbt\squeeze(d(:,m,:)-c);
end
for j=1:p,
    x(:,:,j)=zz*x(:,:,j)*u';
end

% 01/25/03 MJ corrected bug for i==m (sign of c in x determination)
% 01/31/03 MJ added 'real' to qz call
