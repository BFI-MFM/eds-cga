function [G1pi,C,impact,nmat,TT1,TT2,gev,eu, DD, E2, E5, GAMMA, FL_RANK ]=PI_gensys(a0,a1,a2,a3,c,PSI,NX,NETA,FL_RANK,M_,options_)
% System given as
%        a0*E_t[y(t+1])+a1*y(t)=a2*y(t-1)+c+psi*eps(t)
% with z an exogenous variable process.
% Returned system is
%       [s(t)' x(t)' E_t x(t+1)']'=G1pi [s(t-1)' x(t-1)' x(t)]'+C+impact*eps(t),
%  and (a) the matrix nmat satisfying   nmat*E_t z(t)+ E_t x(t+1)=0
%      (b) matrices TT1, TT2  that relate y(t) to these states: y(t)=[TT1 TT2][s(t)' x(t)']'.
% Note that the dimension of the state vector = dim(a0)+NO_FL_EQS
%
% If div is omitted from argument list, a div>1 is calculated.
% eu(1)=1 for existence, eu(2)=1 for uniqueness.  eu(1)=-1 for
% existence only with not-s.c. z; eu=[-2,-2] for coincident zeros.
% Based on
% Christopher A. Sims
% Corrected 10/28/96 by CAS

% Copyright (C) 1996-2009 Christopher Sims
% Copyright (C) 2010-2011 Dynare Team
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


lastwarn('','');
global lq_instruments;
eu=[0;0];C=c;
realsmall=1e-6;
fixdiv=(nargin==6);
n=size(a0,1);
DD=[];E2=[]; E5=0; GAMMA=[];
%
% Find SVD of a0, and create partitions of U, S and V
%
[U0,S0,V0] = svd(a0);

FL_RANK=rank(S0);
U01=U0(1:n,1:FL_RANK);
U02=U0(1:n,FL_RANK+1:n);
V01=V0(1:n,1:FL_RANK);
V02=V0(1:n,FL_RANK+1:n);
S01=S0(1:FL_RANK,1:FL_RANK);

C1=U02'*a1*V01;
C2=U02'*a1*V02;
C3=U02'*a2*V01;
C4=U02'*a2*V02;
C5=U02'*PSI;
Sinv=eye(FL_RANK);
for i=1:FL_RANK
    Sinv(i,i)=1/S01(i,i);
end
F1=Sinv*U01'*a1*V01;
F2=Sinv*U01'*a1*V02;
F3=Sinv*U01'*a2*V01;
F4=Sinv*U01'*a2*V02;
F5=Sinv*U01'*PSI;
singular=0;
warning('', '');
try
    if rcond(C2) < 1e-8
        singular=1;
    else
        warning('off','MATLAB:nearlySingularMatrix');
        warning('off','MATLAB:singularMatrix');
        UAVinv=inv(C2); % i.e. inv(U02'*a1*V02)
        [LastWarningTxt LastWarningID]=lastwarn;
        if any(any(isinf(UAVinv)))==1
            singular=1;
        end
    end
    if singular == 1 || strcmp('MATLAB:nearlySingularMatrix',LastWarningID) == 1 || ...
                 strcmp('MATLAB:illConditionedMatrix',LastWarningID)==1 || ...
                 strcmp('MATLAB:singularMatrix',LastWarningID)==1 
        [C1,C2,C3,C4, C5, F1, F2, F3, F4, F5, M1, M2, UAVinv, FL_RANK, V01, V02] = PI_gensys_singularC(C1,C2,C3,C4, C5, F1, F2, F3, F4, F5, V01, V02, 0);
    end
    warning('on','MATLAB:singularMatrix');
    warning('on','MATLAB:nearlySingularMatrix');
    if (any(any(isinf(UAVinv))) || any(any(isnan(UAVinv)))) 
        if(options_.ACES_solver==1)
            disp('ERROR! saving PI_gensys_data_dump');
            save PI_gensys_data_dump
            error('PI_gensys: Inversion of poss. zero matrix UAVinv=inv(U02''*a1*V02)!');
        else
            warning('PI_gensys: Evading inversion of zero matrix UAVinv=inv(U02''*a1*V02)!');
            eu=[0,0];
            return;
        end
    end
catch
    errmsg=lasterror;
    warning(['error callig PI_gensys_singularC: ' errmsg.message ],'errmsg.identifier');
    %error('errcode',['error callig PI_gensys_singularC: ' errmsg.message ]);
end
%
% Define TT1, TT2
%
TT1=V02;
TT2=V01;

%
%Set up the system matrix for the variables s(t)=V02*Y(t), x(t)=V01*Y(t) and E_t x(t+1)
%                                   and define z(t)'=[s(t)' x(t)]
%
%UAVinv=inv(U02'*a1*V02);
FF=F2; %=Sinv*U01'*a1*V02;
G11=UAVinv*C4;
G12=UAVinv*C3;
G13=-UAVinv*C1;
G31=-FF*G11+F4; % +Sinv*U01'*a2*V02;
G32=-FF*G12+F3; % +Sinv*U01'*a2*V01;
G33=-FF*G13-F1; % -Sinv*U01'*a1*V01;
P1=UAVinv*C5; % *U02'*PSI;
P3=-FF*P1+F5; % + Sinv*U01'*PSI; % This should equal 0
G21=zeros(FL_RANK,(n-FL_RANK));
G22=zeros(FL_RANK,FL_RANK);
G23=eye(FL_RANK);
%H2=zeros(FL_RANK,NX);
num_inst=0;

% New Definitions
Ze11=zeros(NX,NX); 
Ze12=zeros(NX,(n-FL_RANK)); 
Ze134=zeros(NX,FL_RANK);
Ze31=zeros(FL_RANK,NX);

% End of New Definitions

%
% System matrix is called 'G1pi'; Shock matrix is called 'impact'
%

G1pi=[Ze11 Ze12 Ze134 Ze134; P1 G11 G12 G13; Ze31 G21 G22 G23; P3 G31 G32 G33];

impact=[eye(NX,NX); zeros(n+FL_RANK,NX)];

if(options_.ACES_solver==1)
    if isfield(lq_instruments,'names')
        num_inst=size(lq_instruments.names,1);
        if num_inst>0 
            i_var=lq_instruments.inst_var_indices;
            N1=UAVinv*U02'*lq_instruments.B1;
            N3=-FF*N1+Sinv*U01'*lq_instruments.B1;
        else
            error('WARNING: There are no instrumnets for ACES!');
        end
        lq_instruments.N1=N1;
        lq_instruments.N3=N3;
    else
        error('WARNING: There are no instrumnets for ACES!');
    end
    E3=V02*[P1 G11 G12 G13];
    E3= E3+ [zeros(size(V01,1),size(E3,2)-size(V01,2)) V01];
    E2=-E3;
    E5=-V02*N1;
    DD=[zeros(NX,size(N1,2));N1; zeros(FL_RANK,size(N1,2));N3];
    II=eye(num_inst);
    GAMMA=[ E3 -E5 %zeros(size(E3,1),num_inst);
            zeros(num_inst,size(E3,2)), II;
          ];
    eu =[1; 1], nmat=[], gev=[];
    return; % do not check B&K compliancy
end

G0pi=eye(n+FL_RANK+NX);
try
    % In Matlab: [aa bb q z v w] = qz(a,b) s.t. qaz = aa, qbz = bb % 
    % In Octave: [aa bb q z v w] = qz(a,b) s.t. q'az = aa, q'bz=bb %
    % and qzcomplex() extension based on lapack zgges produces same 
    % qz output for Octave as Matlab qz() does for Matlab thus:
    if exist('OCTAVE_VERSION')
        [a b q z]=qzcomplex(G0pi,G1pi);
        q=q';
    else
        [a b q z]=qz(G0pi,G1pi);
    end
catch
    try
        lerror=lasterror;
        disp(['PI_Gensys: ' lerror.message]);
        if 0==strcmp('MATLAB:qz:matrixWithNaNInf',lerror.identifier)
            disp '** Unexpected Error PI_Gensys:qz: ** :';
            button=questdlg('Continue Y/N?','Unexpected Error in qz','No','Yes','Yes'); 
            switch button 
              case 'No' 
                error ('Terminated')
                %case 'Yes'
                
            end
        end
        G1pi=[];impact=[];nmat=[]; gev=[];
        eu=[-2;-2];
        return
    catch
        disp '** Unexpected Error in qz ** :';
        disp lerror.message;
        button=questdlg('Continue Y/N?','Unexpected Error in qz','No','Yes','Yes'); 
        switch button 
          case 'No' 
            error ('Terminated') 
          case 'Yes' 
            G1pi=[];impact=[];nmat=[]; gev=[];
            eu=[-2;-2];
            return
        end
    end
end

if ~fixdiv, div=1.01; end
nunstab=0;
zxz=0;
nn=size(a,1);
for i=1:nn
    % ------------------div calc------------
    if ~fixdiv
        if abs(a(i,i)) > 0
            divhat=abs(b(i,i))/abs(a(i,i));
            % bug detected by Vasco Curdia and Daria Finocchiaro, 2/25/2004  A root of
            % exactly 1.01 and no root between 1 and 1.02, led to div being stuck at 1.01
            % and the 1.01 root being misclassified as stable.  Changing < to <= below fixes this.
            if 1+realsmall<divhat && divhat<=div
                div=.5*(1+divhat);
            end
        end
    end
    % ----------------------------------------
    nunstab=nunstab+(abs(b(i,i))>div*abs(a(i,i)));
    if abs(a(i,i))<realsmall && abs(b(i,i))<realsmall
        zxz=1;
    end
end
div ;
if ~zxz
    [a b q z]=qzdiv(div,a,b,q,z);
end

gev=[diag(a) diag(b)];
if zxz
    disp('Coincident zeros.  Indeterminacy and/or nonexistence.')
    eu=[-2;-2];
    % correction added 7/29/2003.  Otherwise the failure to set output
    % arguments leads to an error message and no output (including eu).
    nmat=[]; %;gev=[]
    return
end
if (FL_RANK ~= nunstab && options_.ACES_solver~=1)
    disp(['Number of unstable variables ' num2str(nunstab)]);
    disp( ['does not match number of expectational equations ' num2str(FL_RANK)]); 
    nmat=[];% gev=[];
    eu=[-2;-2];
    return
end

% New Definitions
z1=z(:,1:n+NX)';
z2=z(:,n+NX+1:n+NX+FL_RANK)';

% New N Matrix by J Pearlman
z12=z2(:,1:n+NX);
z22=z2(:,n+NX+1:n+NX+FL_RANK);
% End of New Definitions

% modified by GP:
nmat=real(inv(z22)*z12);
eu=[1;1];
