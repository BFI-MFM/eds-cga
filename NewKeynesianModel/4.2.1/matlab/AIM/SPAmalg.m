%  [b,rts,ia,nexact,nnumeric,lgroots,aimcode] = ...
%                       SPAmalg(h,neq,nlag,nlead,condn,uprbnd)
%
%  Solve a linear perfect foresight model using the matlab eig
%  function to find the invariant subspace associated with the big
%  roots.  This procedure will fail if the companion matrix is
%  defective and does not have a linearly independent set of
%  eigenvectors associated with the big roots.
% 
%  Input arguments:
% 
%    h         Structural coefficient matrix (neq,neq*(nlag+1+nlead)).
%    neq       Number of equations.
%    nlag      Number of lags.
%    nlead     Number of leads.
%    condn     Zero tolerance used as a condition number test
%              by numeric_shift and reduced_form.
%    uprbnd    Inclusive upper bound for the modulus of roots
%              allowed in the reduced form.
% 
%  Output arguments:
% 
%    b         Reduced form coefficient matrix (neq,neq*nlag).
%    rts       Roots returned by eig.
%    ia        Dimension of companion matrix (number of non-trivial
%              elements in rts).
%    nexact    Number of exact shiftrights.
%    nnumeric  Number of numeric shiftrights.
%    lgroots   Number of roots greater in modulus than uprbnd.
%    aimcode     Return code: see function aimerr.

% Original author: Gary Anderson
% Original file downloaded from:
% http://www.federalreserve.gov/Pubs/oss/oss4/code.html
% Adapted for Dynare by Dynare Team, in order to deal
% with infinite or nan values.
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

function [b,rts,ia,nexact,nnumeric,lgroots,aimcode] = ...
                        SPAmalg(h,neq,nlag,nlead,condn,uprbnd)
b=[];rts=[];ia=[];nexact=[];nnumeric=[];lgroots=[];aimcode=[];
if(nlag<1 || nlead<1) 
    error('Aim_eig: model must have at least one lag and one lead');
end
% Initialization.
nexact=0;nnumeric=0;lgroots=0;iq=0;aimcode=0;b=0;qrows=neq*nlead;qcols=neq*(nlag+nlead);
bcols=neq*nlag;q=zeros(qrows,qcols);rts=zeros(qcols,1);
[h,q,iq,nexact]=SPExact_shift(h,q,iq,qrows,qcols,neq);
if (iq>qrows)
    aimcode = 61;
    return;
end
[h,q,iq,nnumeric]=SPNumeric_shift(h,q,iq,qrows,qcols,neq,condn);
if (iq>qrows)
    aimcode = 62;
    return;
end
[a,ia,js] = SPBuild_a(h,qcols,neq);
if (ia ~= 0)
    if any(any(isnan(a))) || any(any(isinf(a))) 
        disp('A is NAN or INF')
        aimcode=63; 
        return 
    end 
    [w,rts,lgroots,flag_trouble]=SPEigensystem(a,uprbnd,min(length(js),qrows-iq+1));
    if flag_trouble==1; 
        disp('Problem in SPEIG'); 
        aimcode=63;
        return
    end 
    q = SPCopy_w(q,w,js,iq,qrows);
end
test=nexact+nnumeric+lgroots;
if (test > qrows)
    aimcode = 3;
elseif (test < qrows)
    aimcode = 4;
end
if(aimcode==0)
    [nonsing,b]=SPReduced_form(q,qrows,qcols,bcols,neq,condn);
    if ( nonsing && aimcode==0)
        aimcode =  1;
    elseif (~nonsing && aimcode==0)
        aimcode =  5;
    elseif (~nonsing && aimcode==3)
        aimcode = 35;
    elseif (~nonsing && aimcode==4)
        aimcode = 45;
    end
end