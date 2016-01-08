function Z=disc_riccati_fast(F,D,R,H,ch)
% function Z=disc_riccati_fast(F,D,R,H,ch)
% 
% Solves discrete Riccati Equation: 
% Z=FZF' - FZD'inv(DZD'+R)DZF' + H
% Using the Doubling Algorithm 
%
% George Perendia: based on the doubling algorithm for Riccati   
% and the disclyap_fast function provided by Prof. Joe Pearlman 
% V.1 19/5/2006
% V.2 22/10/06
% =================================================================

% Copyright (C) 2006-2011 Dynare Team
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

if nargin == 4 || isempty( ch ) == 1 
    flag_ch = 0; 
else 
    flag_ch = 1; 
end 


% intialisation
tol = 1e-10; % iteration convergence threshold
P0=H; 
X0=F;
if ~any(R) % i.e. ==0
    warning('Dangerously evading inversion of zero matrix!');
    Y0=zeros(size(R));
else
    Y0=D'*inv(R)*D;
end
POYO=P0*Y0;
I=eye(size(POYO));
clear POYO;

% iterate Riccati equation solution
matd=1; 
count=0;
while matd > tol && count < 100
    INVPY=inv(I+P0*Y0);
    P1=X0*INVPY*P0*X0'+ P0; 
    Y1=X0'*Y0*INVPY*X0+ Y0; 
    X1=X0*INVPY*X0; 
    matd=sum( sum(abs( P1 - P0 ))); 
    %    P0=(P1+P1')/2
    P0=P1; 
    X0=X1;
    Y0=Y1;
    count=count+1;
    %    matd;
end 

Z=(P0+P0')/2;
%Z=P0
% check if the convergence took place
if count==100
    matd
    error('Riccati not converged fast enough!');
    %    error.identifier='Riccati not converged!'
    %    error
end
%if count >5 
%    disp('Riccati count= ');
%    count
%end

clear X0 X1 Y0 Y1 P1 I INVPY; 

% Check that X is positive definite 
if flag_ch==1 
    [C,p]=chol(Z); 
    if p ~= 0 
        error('Z is not positive definite')
    end 
end 
