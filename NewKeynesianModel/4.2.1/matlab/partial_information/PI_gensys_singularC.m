function [C1,C2,C3,C4, C5, F1, F2, F3, F4, F5, M1, M2, UAVinv, FL_RANK, V01, V02]=PI_gensys_singularC(C1in, C2in, C3in, C4in, C5in, F1, F2, F3, F4, F5, V01, V02, level)
% [C1,C2,C3,C4, C5, F1, F2, F3, F4, F5, M1, M2, UAVinv,FL_RANK, V01, V02]...
%         =PI_gensys_singularC(C1in, C2in, C3in, C4in, C5in, F1, F2, F3, F4, F5, V01, V02, level)
% 
% Recursive extension for PI_gensys function PCL general DSGE solver
% devised by Prof. Joseph Pearlman 
% developed by George Perendia
% December 2010

% Copyright (C) 1996-2011 Dynare Team
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

level=level+1
if level>100
    error( ' PI_gensys_singularC recurssion exceeeded its maximum of 100 iterations! ');
end
warning('', '');
M1=[];M2=[]; UAVinv=[];
%
% Find SVD of a0, and create partitions of U, S and V
%

[J0,K0,L0] = svd(C2in);
n=size(C2in,1);
K_RANK=rank(K0);
J2=J0(1:n,K_RANK+1:n);

J2C1=J2'*C1in;
M = null(J2C1)';
MJCinv= inv([M;J2C1]);
[sm1, sm2]=size (M');
M1=MJCinv(1:sm1,1:sm2);
M2=MJCinv(1:sm1,1+sm2:end);
FL_RANK=rank(M);

%Define new Cs
C5=[ C5in; J2C1*F5];
C4=[C4in C3in*M2; J2C1*F4 J2C1*F3*M2];
C3=[ C3in*M1; J2C1*F3*M1];
C2=[C2in C1in*M2; J2'*(C4in+C1in*F2) J2'*(C3in+C1in*F1)*M2];
C1=[ C1in*M1; J2'*(C3in+C1in*F1)*M1];
%define new after Cs Fs
% keep this reverse order!!
F5 =M*F5 ;
F4 =[M*F4  M*F3*M2];
F3 = M*F3*M1;
F2 =[M*F2  M*F1*M2];
F1 = M*F1*M1;

V02=[V02 V01*M2];
V01=V01*M1;

warning('', '');
singular=0;
try
    if rcond(C2) < 1e-8
        singular=1;
    else
        UAVinv=inv(C2);
        [LastWarningTxt LastWarningID]=lastwarn;
        if any(any(isinf(UAVinv)))==1
            singular=1;
        end
    end
    % line test is for Octave strncmp('warning: inverse: matrix singular',LastWarningTxt, 33)==1 || ...
    if  singular==1 || strcmp('MATLAB:nearlySingularMatrix',LastWarningID)==1 || ...
                 strcmp('MATLAB:illConditionedMatrix',LastWarningID)==1 || ...
                 strcmp('MATLAB:singularMatrix',LastWarningID)==1
        [C1,C2,C3,C4, C5, F1, F2, F3, F4, F5, M1, M2, UAVinv, FL_RANK, V01, V02] = PI_gensys_singularC(C1,C2,C3,C4, C5, F1, F2, F3, F4, F5, V01, V02, level);
    end
catch
    [errmsg, errcode]=lasterr;
    warning(['error callig PI_gensys_singularC: ' errmsg ],'errcode');
    error('errcode',['error callig PI_gensys_singularC: ' errmsg ]);
end

return;

