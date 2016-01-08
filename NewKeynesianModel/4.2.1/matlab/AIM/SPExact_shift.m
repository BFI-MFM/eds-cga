% [h,q,iq,nexact] = exact_shift(h,q,iq,qrows,qcols,neq)
%
% Compute the exact shiftrights and store them in q.

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

function [h,q,iq,nexact] = SPExact_shift(h,q,iq,qrows,qcols,neq)

%hs=SPSparse(h);
hs=sparse(h);
nexact = 0;
left   = 1:qcols;
right  = qcols+1:qcols+neq;
zerorows = find( sum(abs( hs(:,right)' ))==0 );

while( any(zerorows) && iq <= qrows )
   nz = length(zerorows);
   q(iq+1:iq+nz,:) = hs(zerorows,left);
   hs(zerorows,:)   = SPShiftright(hs(zerorows,:),neq);
   iq     = iq + nz;
   nexact = nexact + nz;
   zerorows = find( sum(abs( hs(:,right)' ))==0 );
end
h=full(hs);

