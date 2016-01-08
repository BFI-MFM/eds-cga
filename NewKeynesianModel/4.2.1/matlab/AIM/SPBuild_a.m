%  [a,ia,js] = SPBuild_a(h,qcols,neq)
%
%  Build the companion matrix, deleting inessential lags.
%  Solve for x_{t+nlead} in terms of x_{t+nlag},...,x_{t+nlead-1}.

% Original author: Gary Anderson
% Original file downloaded from:
% http://www.federalreserve.gov/Pubs/oss/oss4/code.html
% Adapted for Dynare by Dynare Team.
%
% This code in the public domain and may be used freely.
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

function [a,ia,js] = SPBuild_a(h,qcols,neq)

left  = 1:qcols;
right = qcols+1:qcols+neq;
%hs=SPSparse(h);
hs=sparse(h);

hs(:,left) = -hs(:,right)\hs(:,left);

%  Build the big transition matrix.

a = zeros(qcols,qcols);

if(qcols > neq)
  eyerows = 1:qcols-neq;
  eyecols = neq+1:qcols;
  a(eyerows,eyecols) = eye(qcols-neq);
end
hrows      = qcols-neq+1:qcols;
a(hrows,:) = hs(:,left);

%  Delete inessential lags and build index array js.  js indexes the
%  columns in the big transition matrix that correspond to the
%  essential lags in the model.  They are the columns of q that will
%  get the unstable left eigenvectors. 

js       = 1:qcols;
zerocols = sum(abs(a)) == 0;
while( any(zerocols) )
  a(:,zerocols) = [];
  a(zerocols,:) = [];
  js(zerocols)  = [];
  zerocols = sum(abs(a)) == 0;
end
ia = length(js);
