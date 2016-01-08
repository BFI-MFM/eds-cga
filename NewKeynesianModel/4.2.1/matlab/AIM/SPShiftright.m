function [y] = SPShiftright(x,n)

% [y] = shiftright(x,n)
%
%  Shift the rows of x to the right by n columns, leaving zeros in the
%  first n columns. 

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

[rows,cols] = size(x);

left  = 1:cols-n;
right = n+1:cols;

y = zeros(rows,cols);
y(:,right) = x(:,left);

return
