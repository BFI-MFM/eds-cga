function [dr,aimcode,rts]=dynAIMsolver1(jacobia_,M_,dr)
% function [dr,aimcode]=dynAIMsolver1(jacobia_,M_,dr)
% Maps Dynare jacobia to AIM 1st order model solver designed and developed by Gary ANderson 
% and derives the solution for gy=dr.hgx and gu=dr.hgu from the AIM outputs
% AIM System is given as a sum: 
% i.e. for i=-$...+&   SUM(Hi*xt+i)= £*zt, t = 0, . . . ,?
% and its input as single array of matrices: [H-$...  Hi ... H+&]
% and its solution as xt=SUM( Bi*xt+i) + @*£*zt for i=-$...-1 
% with the output in form bb=[B-$...  Bi ... B-1] and @=inv(Ho+H1*B-1) 
% Dynare jacobian = [fy'-$...  fy'i ... fy'+&  fu'] 
% where [fy'-$...  fy'i ... fy'+&]=[H-$...  Hi ... H+&] and fu'= £
%
% INPUTS
%   jacobia_   [matrix]           1st order derivative of the model 
%   dr         [matlab structure] Decision rules for stochastic simulations.
%   M_         [matlab structure] Definition of the model.           
%    
% OUTPUTS
%   dr         [matlab structure] Decision rules for stochastic simulations.
%   aimcode    [integer]          1: the model defines variables uniquely
%   aimcode is resolved in AIMerr as
%      (c==1)  e='Aim: unique solution.';
%      (c==2)  e='Aim: roots not correctly computed by real_schur.';
%      (c==3)  e='Aim: too many big roots.';
%      (c==35) e='Aim: too many big roots, and q(:,right) is singular.';
%      (c==4)  e='Aim: too few big roots.';
%      (c==45) e='Aim: too few big roots, and q(:,right) is singular.';
%      (c==5)  e='Aim: q(:,right) is singular.';
%      (c==61) e='Aim: too many exact shiftrights.';
%      (c==62) e='Aim: too many numeric shiftrights.';
%      else    e='Aimerr: return code not properly specified';
%    
% SPECIAL REQUIREMENTS
% Dynare use: 
%       1) the lognormal block in DR1 is being invoked for some models and changing
%       values of ghx and ghy. We need to return the AIM output
%       values before that block and run the block with the current returned values
%       of gy (i.e. dr.ghx) and gu (dr.ghu) if it is needed even when the AIM is used  
%       (it does not depend on mjdgges output).
%       
%       2) passing in aa={Q'|1}*jacobia_ can produce ~ one order closer  
%       results to the Dynare solutiion then when if plain jacobia_ is passed, 
%       i.e. diff < e-14 for aa and diff < *e-13 for jacobia_ if Q' is used.  
%
% GP July 2008 

% Copyright (C) 2008-2009 Dynare Team
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

aimcode=-1;
neq= size(jacobia_,1); % no of equations
lags=M_.maximum_endo_lag; % no of lags and leads
leads=M_.maximum_endo_lead;
klen=(leads+lags+1);  % total lenght
theAIM_H=zeros(neq, neq*klen); % alloc space
lli=M_.lead_lag_incidence';  
% "sparse" the compact jacobia into AIM H aray of matrices 
% without exogenous shoc
theAIM_H(:,find(lli(:)))=jacobia_(:,nonzeros(lli(:)));
condn  = 1.e-10;%SPAmalg uses this in zero tests
uprbnd = 1 + 1.e-6;%allow unit roots
                   % forward only models - AIM must have at least 1 lead and 1 lag.
if lags ==0
    theAIM_H =[zeros(neq) theAIM_H];
    lags=1;
    klen=(leads+lags+1);
end
% backward looking only models
if leads ==0
    theAIM_H =[theAIM_H zeros(neq)];
    leads=1;
    klen=(leads+lags+1);
end
%disp('gensysToAMA:running ama');
try % try to run AIM
    [bb,rts,ia,nexact,nnumeric,lgroots,aimcode] =...
        SPAmalg(theAIM_H,neq, lags,leads,condn,uprbnd);
catch
    err = lasterror;
    disp(['Dynare AIM Solver error:' sprintf('%s; ID:%s',err.message, err.identifier)]);
    rethrow(lasterror);
end
if aimcode==1 %if OK
    col_order=[];
    for i =1:lags
        col_order=[((i-1)*neq)+dr.order_var' col_order];
    end
    bb_ord= bb(dr.order_var,col_order); % derive ordered gy
    
    % variables are present in the state space at the lag at which they
    % appear and at all smaller lags. The are ordered from smaller to
    % higher lag (reversed order of M_.lead_lag_incidence rows for lagged
    % variables)
    i_lagged_vars = flipud(cumsum(M_.lead_lag_incidence(1:M_.maximum_lag,dr.order_var),1))';

    dr.ghx= bb_ord(:,find(i_lagged_vars(:))); % get columns reported in
                                              % Dynare solution
    if M_.exo_nbr % if there are exogenous shocks then derive gu for the shocks:
                  %   get H0 and H+1=HM
                  %    theH0= theAIM_H (:,M_.maximum_endo_lag*neq+1: (M_.maximum_endo_lag+1)*neq);
                  %theH0= theAIM_H (:,lags*neq+1: (lags+1)*neq);
                  %    theHP= theAIM_H (:,(M_.maximum_endo_lag+1)*neq+1: (M_.maximum_endo_lag+2)*neq);
                  %theHP= theAIM_H (:,(lags+1)*neq+1: (lags+2)*neq);
        theAIM_Psi= - jacobia_(:, size(nonzeros(lli(:)))+1:end);% 
                                                                %? = inv(H0 + H1B1)
                                                                %phi= (theH0+theHP*sparse(bb(:,(lags-1)*neq+1:end)))\eye( neq);
                                                                %AIM_ghu=phi*theAIM_Psi;
                                                                %dr.ghu =AIM_ghu(dr.order_var,:); % order gu
                                                                % Using AIM SPObstruct
        scof = SPObstruct(theAIM_H,bb,neq,lags,leads);
        scof1= scof(:,(lags)*neq+1:end);
        scof1= scof1(:,dr.order_var);
        dr.ghu =scof1\theAIM_Psi;
    else
        dr.ghu = [];
    end
else
    err=SPAimerr(aimcode);
    %warning('Error in AIM: aimcode=%d, erro=%s', aimcode, err);;
    disp(['Error in AIM: aimcode=' sprintf('%d : %s',aimcode, err)]);
    if aimcode < 1 || aimcode > 5  % too big exception, use mjdgges
        error('Error in AIM: aimcode=%d ; %s', aimcode, err);
    end
    %    if aimcode > 5 
    %        disp(['Error in AIM: aimcode=' sprintf('%d : %s',aimcode, err)]);
    %        aimcode=5;
    %    end  
end
