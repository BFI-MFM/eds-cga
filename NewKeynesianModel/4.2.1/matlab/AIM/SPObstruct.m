% scof = SPObstruct(cof,cofb,neq,nlag,nlead)
%
% Construct the coefficients in the observable structure.
%    
%   Input arguments:
%
%            cof    structural coefficients
%            cofb   reduced form
%            neq    number of equations
%            nlag   number of lags
%            nlead  number of leads
%
%   Output arguments:
%
%            scof  observable structure coefficients

% Original author: Gary Anderson
% Original file downloaded from:
% http://www.federalreserve.gov/Pubs/oss/oss4/code.html
% Adapted for Dynare by Dynare Team.
%
% This code is in the public domain and may be used freely.
% However the authors would appreciate acknowledgement of the source by
% citation of any of the following papers:
%
% Anderson, G. and Moore, G.
% "A Linear Algebraic Procedure for Solving Linear Perfect Foresight
% Models."
% Economics Letters, 17, 1985.
%
% Anderson, G.
% "Solving Linear Rational Expectations Models: A Horse Race"
% Computational Economics, 2008, vol. 31, issue 2, pages 95-113
%
% Anderson, G.
% "A Reliable and Computationally Efficient Algorithm for Imposing the
% Saddle Point Property in Dynamic Models"
% Journal of Economic Dynamics and Control, Forthcoming

function scof = SPObstruct(cof,cofb,neq,nlag,nlead)


% Append the negative identity to cofb

cofb = [cofb, -eye(neq)];
scof = zeros(neq,neq*(nlag+1));
q    = zeros(neq*nlead, neq*(nlag+nlead));
[rc,cc] = size(cofb);
qs=sparse(q);
qs(1:rc,1:cc) = sparse(cofb);
qcols = neq*(nlag+nlead);

if( nlead > 1 ) 
   for i = 1:nlead-1
      rows = i*neq + (1:neq);
      qs(rows,:) = SPShiftright( qs((rows-neq),:), neq );
   end
end

l = (1: neq*nlag);
r = (neq*nlag+1: neq*(nlag+nlead));

 qs(:,l) = - qs(:,r) \ qs(:,l);

minus =              1:       neq*(nlag+1);
plus  = neq*(nlag+1)+1: neq*(nlag+1+nlead);

cofs=sparse(cof);
scof(:,neq+1:neq*(nlag+1)) = cofs(:,plus)*qs(:,l);
scof = scof + cofs(:,minus);
