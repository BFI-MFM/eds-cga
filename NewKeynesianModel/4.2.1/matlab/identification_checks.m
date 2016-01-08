function [McoH, McoJ, McoGP, PcoH, PcoJ, PcoGP, condH, condJ, condGP, eH, eJ, eGP, ind01, ind02, indnoH, indnoJ, ixnoH, ixnoJ] = identification_checks(H, JJ, gp)
% function [McoH, McoJ, McoGP, PcoH, PcoJ, PcoGP, condH, condJ, condGP, eH,
%          eJ, eGP, ind01, ind02, indnoH, indnoJ, ixnoH, ixnoJ] = identification_checks(H, JJ, gp)
% checks for identification
%
% INPUTS
%    o H                [matrix] [(entries in st.sp. model solutio) x nparams] 
%                                derivatives of model solution w.r.t. parameters and shocks
%    o JJ               [matrix] [moments x nparams] 
%                                 derivatives of moments w.r.t. parameters and shocks
%    o gp               [matrix] [jacobian_entries x nparams] 
%                                derivatives of jacobian (i.e. LRE model) w.r.t. parameters and shocks
%    
% OUTPUTS
%    o McoH             [array] multicollinearity coefficients in the model solution
%    o McoJ             [array] multicollinearity coefficients in the moments
%    o McoGP            [array] multicollinearity coefficients in the LRE model
%    o PcoH             [matrix] pairwise correlations in the model solution
%    o PcoJ             [matrix] pairwise correlations in the moments
%    o PcoGP            [matrix] pairwise correlations in the LRE model
%    o condH            condition number of H
%    o condJ            condition number of JJ
%    o condGP           condition number of gp
%    o eH               eigevectors of H
%    o eJ               eigevectors of JJ
%    o eGP              eigevectors of gp
%    o ind01            [array] binary indicator for zero columns of H
%    o ind02            [array] binary indicator for zero columns of JJ
%    o indnoH           [matrix] index of non-identified params in H
%    o indnoJ           [matrix] index of non-identified params in JJ
%    o ixnoH            number of rows in ind01
%    o ixnoJ            number of rows in ind02 
%    
% SPECIAL REQUIREMENTS
%    None

% Copyright (C) 2008-2011 Dynare Team
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

% My suggestion is to have the following steps for identification check in
% dynare:

% 1. check rank of H, JJ, gp at theta
npar = size(H,2);
npar0 = size(gp,2);  % shocks do not enter jacobian
indnoH = zeros(1,npar);
indnoJ = zeros(1,npar);
indnoLRE = zeros(1,npar0);

% H matrix
ind1 = find(vnorm(H)>=eps); % take non-zero columns
H1 = H(:,ind1);
[eu,e2,e1] = svd( H1, 0 );
eH = zeros(npar,npar);
% eH(ind1,:) = e1;
eH(ind1,length(find(vnorm(H)<eps))+1:end) = e1; % non-zero eigenvectors
eH(find(vnorm(H)<eps),1:length(find(vnorm(H)<eps)))=eye(length(find(vnorm(H)<eps)));
condH = cond(H1);
rankH = rank(H);
rankHH = rank(H'*H);

ind2 = find(vnorm(JJ)>=eps); % take non-zero columns
JJ1 = JJ(:,ind2);
[eu,ee2,ee1] = svd( JJ1, 0 );
eJ = zeros(npar,npar);
eJ(ind2,length(find(vnorm(JJ)<eps))+1:end) = ee1; % non-zero eigenvectors
eJ(find(vnorm(JJ)<eps),1:length(find(vnorm(JJ)<eps)))=eye(length(find(vnorm(JJ)<eps)));
condJ = cond(JJ1);
rankJJ = rank(JJ'*JJ);
rankJ = rank(JJ);

ind3 = find(vnorm(gp)>=eps); % take non-zero columns
gp1 = gp(:,ind3);
covgp = gp1'*gp1;
sdgp = sqrt(diag(covgp));
sdgp = sdgp*sdgp';
[eu,ex2,ex1] = svd(gp1, 0 );
eGP = zeros(npar0,npar0);
eGP(ind3,length(find(vnorm(gp)<eps))+1:end) = ex1; % non-zero eigenvectors
eGP(find(vnorm(gp)<eps),1:length(find(vnorm(gp)<eps)))=eye(length(find(vnorm(gp)<eps)));
% condJ = cond(JJ1'*JJ1);
condGP = cond(gp1);


ind01 = zeros(npar,1);
ind02 = zeros(npar,1);
ind01(ind1) = 1;
ind02(ind2) = 1;

% find near linear dependence problems:
McoH = NaN(npar,1);
McoJ = NaN(npar,1);
McoGP = NaN(npar0,1);
for ii = 1:size(H1,2);
    McoH(ind1(ii),:) = [cosn([H1(:,ii),H1(:,find([1:1:size(H1,2)]~=ii))])];
end
for ii = 1:size(JJ1,2);
    McoJ(ind2(ii),:) = [cosn([JJ1(:,ii),JJ1(:,find([1:1:size(JJ1,2)]~=ii))])];
end
for ii = 1:size(gp1,2);
    McoGP(ind3(ii),:) = [cosn([gp1(:,ii),gp1(:,find([1:1:size(gp1,2)]~=ii))])];
end


ixno = 0;
if rankH<npar | rankHH<npar | min(1-McoH)<1.e-10
    %         - find out which parameters are involved,
    %   using the vnorm and the svd of H computed before;
    %   disp('Some parameters are NOT identified in the model: H rank deficient')
    %   disp(' ')
    if length(ind1)<npar,
        % parameters with zero column in H
        ixno = ixno + 1;
        indnoH(ixno,:) = (~ismember([1:npar],ind1));
    end
    e0 = [rankHH+1:length(ind1)];
    for j=1:length(e0),
        % linearely dependent parameters in H
        ixno = ixno + 1;
        indnoH(ixno,ind1) = (abs(e1(:,e0(j))) > 1.e-6 )';
    end
else % rank(H)==length(theta), go to 2
     % 2. check rank of J
     %   disp('All parameters are identified at theta in the model (rank of H)')
     %   disp(' ')
end
ixnoH=ixno;

ixno = 0;
if rankJ<npar | rankJJ<npar | min(1-McoJ)<1.e-10
    %         - find out which parameters are involved
    %   disp('Some parameters are NOT identified by the moments included in J')
    %   disp(' ')
    if length(ind2)<npar,
        % parameters with zero column in JJ
        ixno = ixno + 1;
        indnoJ(ixno,:) = (~ismember([1:npar],ind2));
    end
    ee0 = [rankJJ+1:length(ind2)];
    for j=1:length(ee0),
        % linearely dependent parameters in JJ
        ixno = ixno + 1;
        indnoJ(ixno,ind2) = (abs(ee1(:,ee0(j))) > 1.e-6)';
    end
else  %rank(J)==length(theta) =>
      %         disp('All parameters are identified at theta by the moments included in J')
end
ixnoJ=ixno;

% here there is no exact linear dependence, but there are several
%     near-dependencies, mostly due to strong pairwise colliniearities, which can
%     be checked using

PcoH = NaN(npar,npar);
PcoJ = NaN(npar,npar);
PcoGP = NaN(npar0,npar0);
for ii = 1:size(H1,2);
    PcoH(ind1(ii),ind1(ii)) = 1;
    for jj = ii+1:size(H1,2);
        PcoH(ind1(ii),ind1(jj)) = [cosn([H1(:,ii),H1(:,jj)])];
        PcoH(ind1(jj),ind1(ii)) = PcoH(ind1(ii),ind1(jj));
    end
end

for ii = 1:size(JJ1,2);
    PcoJ(ind2(ii),ind2(ii)) = 1;
    for jj = ii+1:size(JJ1,2);
        PcoJ(ind2(ii),ind2(jj)) = [cosn([JJ1(:,ii),JJ1(:,jj)])];
        PcoJ(ind2(jj),ind2(ii)) = PcoJ(ind2(ii),ind2(jj));
    end
end

for ii = 1:size(gp1,2);
    PcoGP(ind3(ii),ind3(ii)) = 1;
    for jj = ii+1:size(gp1,2);
        PcoGP(ind3(ii),ind3(jj)) = [cosn([gp1(:,ii),gp1(:,jj)])];
        PcoGP(ind3(jj),ind3(ii)) = PcoGP(ind3(ii),ind3(jj));
    end
end





