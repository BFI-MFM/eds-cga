% [nonsing,b] = SPReduced_form(q,qrows,qcols,bcols,neq,b,condn);
%
% Compute reduced-form coefficient matrix, b.

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

function [nonsing,b] = SPReduced_form(q,qrows,qcols,bcols,neq,condn);
b=[];
%qs=SPSparse(q);
qs=sparse(q);
left = 1:qcols-qrows;
right = qcols-qrows+1:qcols;
nonsing = rcond(full(qs(:,right))) > condn;
if(nonsing)
    qs(:,left) = -qs(:,right)\qs(:,left);
    b = qs(1:neq,1:bcols);
    b = full(b);
else  %rescale by dividing row by maximal qr element
    %'inverse condition number small, rescaling '
    themax=max(abs(qs(:,right)),[],2);
    oneover = diag(1 ./ themax);
    nonsing = rcond(full(oneover *qs(:,right))) > condn;
    if(nonsing)
        qs(:,left) = -(oneover*qs(:,right))\(oneover*(qs(:,left)));
        b = qs(1:neq,1:bcols);
        b = full(b);
    end
end

