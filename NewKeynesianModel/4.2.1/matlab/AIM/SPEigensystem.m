%  [w,rts,lgroots] = SPEigensystem(a,uprbnd)
%
%  Compute the roots and the left eigenvectors of the companion
%  matrix, sort the roots from large-to-small, and sort the
%  eigenvectors conformably.  Map the eigenvectors into the real
%  domain. Count the roots bigger than uprbnd.

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

function [w,rts,lgroots,flag_trouble] = SPEigensystem(a,uprbnd,rowsLeft) 
opts.disp=0; 
% next block is commented out because eigs() intermitently returns different rts
%try
%    [w,d]   = eigs(a',rowsLeft,'LM',opts);
%    rts     = diag(d);
%    mag     = abs(rts);
%    [mag,k] = sort(-mag);
%    rts     = rts(k);
%catch
    %disp('Catch in SPE');
    %pause(0.5);
    %aStr=datestr(clock);
    %eval(['save ' regexprep(aStr,' ','')  ' a']);
    try
        [w,d]=eig(a');
    catch
        lasterr
        w=[];rts=[];lgroots=[];
        flag_trouble=1;
        return
    end
    rts     = diag(d);
    mag     = abs(rts);
    [mag,k] = sort(-mag);
    rts     = rts(k);
%end
flag_trouble=0; 

%ws=SPSparse(w);
ws=sparse(w);
ws       = ws(:,k);

%  Given a complex conjugate pair of vectors W = [w1,w2], there is a
%  nonsingular matrix D such that W*D = real(W) + imag(W).  That is to
%  say, W and real(W)+imag(W) span the same subspace, which is all
%  that aim cares about. 

ws = real(ws) + imag(ws);

lgroots = sum(abs(rts) > uprbnd);

w=full(ws);