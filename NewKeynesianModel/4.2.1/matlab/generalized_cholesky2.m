function AA = generalized_cholesky2(A)
%function AA = generalized_cholesky2(A)
%
% This procedure produces:
%
% y = chol(A+E), where E is a diagonal matrix with each element as small
% as possible, and A and E are the same size.  E diagonal values are 
% constrained by iteravely updated Gerschgorin bounds.  
%
% REFERENCES:
%
% Jeff Gill and Gary King. 1998. "`Hessian not Invertable.' Help!"
% manuscript in progress, Harvard University.
%
% Robert B. Schnabel and Elizabeth Eskow. 1990. "A New Modified Cholesky
% Factorization," SIAM Journal of Scientific Statistical Computating,
% 11, 6: 1136-58.

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

n = size(A,1);
L = zeros(n,n);
deltaprev = 0;
gamm = max(abs(diag(A))); 
tau = eps^(1/3);

if  min(eig(A)) > 0;
    tau = -1000000;
end;

norm_A = max(sum(abs(A))');
gamm = max(abs(diag(A))); 
delta = max([eps*norm_A;eps]);
Pprod = eye(n);

if n > 2; 
    for k = 1,n-2;
        if min((diag(A(k+1:n,k+1:n))' - A(k,k+1:n).^2/A(k,k))') < tau*gamm ...
                & min(eig(A((k+1):n,(k+1):n))) < 0;
            [tmp,dmax] = max(diag(A(k:n,k:n)));
            if A(k+dmax-1,k+dmax-1) > A(k,k);
                P = eye(n);
                Ptemp = P(k,:); 
                P(k,:) = P(k+dmax-1,:); 
                P(k+dmax-1,:) = Ptemp;
                A = P*A*P;
                L = P*L*P;
                Pprod = P*Pprod;
            end;
            g = zeros(n-(k-1),1);
            for i=k:n;  
                if i == 1;
                    sum1 = 0;
                else;
                    sum1 = sum(abs(A(i,k:(i-1)))');
                end;
                if i == n;
                    sum2 = 0;
                else;
                    sum2 = sum(abs(A((i+1):n,i)));
                end; 
                g(i-k+1) = A(i,i) - sum1 - sum2;
            end; 
            [tmp,gmax] = max(g);
            if gmax ~= k;
                P = eye(n);
                Ptemp  = P(k,:); 
                P(k,:) = P(k+dmax-1,:); 
                P(k+dmax-1,:) = Ptemp;
                A = P*A*P;
                L = P*L*P;
                Pprod = P*Pprod;
            end; 
            normj = sum(abs(A(k+1:n,k)));
            delta = max([0;deltaprev;-A(k,k)+normj;-A(k,k)+tau*gamm]);
            if delta > 0;
                A(k,k) = A(k,k) + delta;
                deltaprev = delta;
            end;
        end; 
        A(k,k) = sqrt(A(k,k));
        L(k,k) = A(k,k); 
        for i=k+1:n; 
            if L(k,k) > eps; 
                A(i,k) = A(i,k)/L(k,k); 
            end;
            L(i,k) = A(i,k);
            A(i,k+1:i) = A(i,k+1:i) - L(i,k)*L(k+1:i,k)';
            if A(i,i) < 0; 
                A(i,i) = 0; 
            end;
        end;
    end;
end;
A(n-1,n) = A(n,n-1);
eigvals  = eig(A(n-1:n,n-1:n));
dlist    = [ 0 ; deltaprev ;...
             -min(eigvals)+tau*max((inv(1-tau)*max(eigvals)-min(eigvals))|gamm) ]; 
if dlist(1) > dlist(2); 
    delta = dlist(1);   
else;
    delta = dlist(2);
end;
if delta < dlist(3);
    delta = dlist(3);
end;
if delta > 0;
    A(n-1,n-1) = A(n-1,n-1) + delta;
    A(n,n) = A(n,n) + delta;
    deltaprev = delta;
end;
A(n-1,n-1) = sqrt(A(n-1,n-1));
L(n-1,n-1) = A(n-1,n-1);
A(n,n-1) = A(n,n-1)/L(n-1,n-1);
L(n,n-1) = A(n,n-1);
A(n,n) = sqrt(A(n,n) - L(n,(n-1))^2);
L(n,n) = A(n,n);
AA = Pprod'*L'*Pprod';