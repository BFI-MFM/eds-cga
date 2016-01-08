% q = SPCopy_w(q,w,js,iq,qrows)
%
%  Copy the eigenvectors corresponding to the largest roots into the
%  remaining empty rows and columns js of q 

% Author: Gary Anderson
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

function  q = SPCopy_w(q,w,js,iq,qrows)

if(iq < qrows)
   lastrows = iq+1:qrows;
   wrows    = 1:length(lastrows);
   q(lastrows,js) = w(:,wrows)';
end
