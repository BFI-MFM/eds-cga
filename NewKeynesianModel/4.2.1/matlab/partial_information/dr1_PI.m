function [dr,info,M_,options_,oo_] = dr1_PI(dr,task,M_,options_,oo_)
% function [dr,info,M_,options_,oo_] = dr1_PI(dr,task,M_,options_,oo_)
% Computes the reduced form solution of a rational expectation model first
% order
% approximation of the Partial Information stochastic model solver around the deterministic steady state).
% Prepares System as
%        A0*E_t[y(t+1])+A1*y(t)=A2*y(t-1)+c+psi*eps(t)
% with z an exogenous variable process.
% and calls PI_Gensys.m solver 
% based on Pearlman et al 1986 paper and derived from
% C.Sims' gensys linear solver.
% to return solution in format
%       [s(t)' x(t)' E_t x(t+1)']'=G1pi [s(t-1)' x(t-1)' x(t)]'+C+impact*eps(t),
%
% INPUTS
%   dr         [matlab structure] Decision rules for stochastic simulations.
%   task       [integer]          if task = 0 then dr1 computes decision rules.
%                                 if task = 1 then dr1 computes eigenvalues.
%   M_         [matlab structure] Definition of the model.           
%   options_   [matlab structure] Global options.
%   oo_        [matlab structure] Results 
%    
% OUTPUTS
%   dr         [matlab structure] Decision rules for stochastic simulations.
%   info       [integer]          info=1: the model doesn't define current variables uniquely
%                                 info=2: problem in mjdgges.dll info(2) contains error code. 
%                                 info=3: BK order condition not satisfied info(2) contains "distance"
%                                         absence of stable trajectory.
%                                 info=4: BK order condition not satisfied info(2) contains "distance"
%                                         indeterminacy.
%                                 info=5: BK rank condition not satisfied.
%   M_         [matlab structure]            
%   options_   [matlab structure]
%   oo_        [matlab structure]
%  
% ALGORITHM
%   ...
%    
% SPECIAL REQUIREMENTS
%   none.
%  

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

global lq_instruments;
info = 0;

options_ = set_default_option(options_,'loglinear',0);
options_ = set_default_option(options_,'noprint',0);
options_ = set_default_option(options_,'olr',0);
options_ = set_default_option(options_,'olr_beta',1);
options_ = set_default_option(options_,'qz_criterium',1.000001);

xlen = M_.maximum_endo_lead + M_.maximum_endo_lag + 1;

if (options_.aim_solver == 1)
    options_.aim_solver == 0;
    warning('You can not use AIM with Part Info solver. AIM ignored'); 
end
if (options_.order > 1) 
    warning('You can not use order higher than 1 with Part Info solver. Order 1 assumed');
    options_.order =1;
end

% expanding system for Optimal Linear Regulator
if options_.ramsey_policy && options_.ACES_solver == 0
    if isfield(M_,'orig_model')
        orig_model = M_.orig_model;
        M_.endo_nbr = orig_model.endo_nbr;
        M_.endo_names = orig_model.endo_names;
        M_.lead_lag_incidence = orig_model.lead_lag_incidence;
        M_.maximum_lead = orig_model.maximum_lead;
        M_.maximum_endo_lead = orig_model.maximum_endo_lead;
        M_.maximum_lag = orig_model.maximum_lag;
        M_.maximum_endo_lag = orig_model.maximum_endo_lag;
    end
    old_solve_algo = options_.solve_algo;
    %  options_.solve_algo = 1;
    oo_.steady_state = dynare_solve('ramsey_static',oo_.steady_state,0,M_,options_,oo_,it_);
    options_.solve_algo = old_solve_algo;
    [junk,junk,multbar] = ramsey_static(oo_.steady_state,M_,options_,oo_,it_);
    [jacobia_,M_] = ramsey_dynamic(oo_.steady_state,multbar,M_,options_,oo_,it_);
    klen = M_.maximum_lag + M_.maximum_lead + 1;
    dr.ys = [oo_.steady_state;zeros(M_.exo_nbr,1);multbar];
else
    klen = M_.maximum_lag + M_.maximum_lead + 1;
    iyv = M_.lead_lag_incidence';
    iyv = iyv(:);
    iyr0 = find(iyv) ;
    it_ = M_.maximum_lag + 1 ;
    
    if M_.exo_nbr == 0
        oo_.exo_steady_state = [] ;
    end
    

    if options_.ACES_solver == 1
        sim_ruleids=[];
        tct_ruleids=[];
        if  size(M_.equations_tags,1)>0  % there are tagged equations, check if they are aceslq rules
            for teq=1:size(M_.equations_tags,1)
                if strcmp(M_.equations_tags(teq,3),'aceslq_sim_rule')
                    sim_ruleids=[sim_ruleids cell2mat(M_.equations_tags(teq,1))]
                end
                if strcmp(M_.equations_tags(teq,3),'aceslq_tct_rule')
                    tct_ruleids=[tct_ruleids cell2mat(M_.equations_tags(teq,1))]
                end
            end
        end
        lq_instruments.sim_ruleids=sim_ruleids;
        lq_instruments.tct_ruleids=tct_ruleids;
        %if isfield(lq_instruments,'xsopt_SS') %% changed by BY
        [junk, lq_instruments.xsopt_SS,lq_instruments.lmopt_SS,s2,check] = opt_steady_get;%% changed by BY
        [qc, DYN_Q] = QPsolve(lq_instruments, s2, check); %% added by BY
        z = repmat(lq_instruments.xsopt_SS,1,klen);
    else
        z = repmat(dr.ys,1,klen);
    end
    z = z(iyr0) ;
    [junk,jacobia_] = feval([M_.fname '_dynamic'],z,[oo_.exo_simul ...
                        oo_.exo_det_simul], M_.params, dr.ys, it_);

    if options_.ACES_solver==1 && (length(sim_ruleids)>0 || length(tct_ruleids)>0 )
        if length(sim_ruleids)>0
            sim_rule=jacobia_(sim_ruleids,:);
            % uses the subdirectory - BY
            save ([ACES_DirectoryName,'/',M_.fname '_ACESLQ_sim_rule.txt'], 'sim_rule', '-ascii', '-double', '-tabs');
        end
        if length(tct_ruleids)>0
            tct_rule=jacobia_(tct_ruleids,:);
            % uses the subdirectory - BY
            save ([ACES_DirectoryName,'/',M_.fname '_ACESLQ_tct_rule.txt'], 'tct_rule', '-ascii', '-double', '-tabs');
        end
        aces_ruleids=union(tct_ruleids,sim_ruleids);
        j_size=size(jacobia_,1);
        j_rows=1:j_size;
        j_rows = setxor(j_rows,aces_ruleids);
        jacobia_=jacobia_(j_rows ,:);
    end

end

if options_.debug
    save([M_.fname '_debug.mat'],'jacobia_')
end

dr=set_state_space(dr,M_);
kstate = dr.kstate;
kad = dr.kad;
kae = dr.kae;
nstatic = dr.nstatic;
nfwrd = dr.nfwrd;
npred = dr.npred;
nboth = dr.nboth;
order_var = dr.order_var;
nd = size(kstate,1);
nz = nnz(M_.lead_lag_incidence);

sdyn = M_.endo_nbr - nstatic;

k0 = M_.lead_lag_incidence(M_.maximum_endo_lag+1,order_var);
k1 = M_.lead_lag_incidence(find([1:klen] ~= M_.maximum_endo_lag+1),:);

if (options_.aim_solver == 1) 
    error('Anderson and Moore AIM solver is not compatible with Partial Information models');
end % end if useAIM and...

    % If required, try PCL86 solver, that is, if not the check being
    % performed only and if it is 1st order 
    % create sparse, extended jacobia AA:
    nendo=M_.endo_nbr; % = size(aa,1)
    

    if(options_.ACES_solver==1)
        %if ~isfield(lq_instruments,'names')
        if isfield(options_,'instruments')
            lq_instruments.names=options_.instruments;
        end
        %end
        if isfield(lq_instruments,'names')
            num_inst=size(lq_instruments.names,1);
            if ~isfield(lq_instruments,'inst_var_indices') && num_inst>0
                for i=1:num_inst
                    i_tmp = strmatch(deblank(lq_instruments.names(i,:)),M_.endo_names,'exact');
                    if isempty(i_tmp)
                        error (['One of the specified instrument variables does not exist']) ;
                    else
                        i_var(i) = i_tmp;
                    end
                end
                lq_instruments.inst_var_indices=i_var;
            elseif size(lq_instruments.inst_var_indices)>0
                i_var=lq_instruments.inst_var_indices;
                if ~num_inst
                    num_inst=size(lq_instruments.inst_var_indices);
                end
            else
                i_var=[];
                num_inst=0;
            end
            if size(i_var,2)>0 && size(i_var,2)==num_inst
                m_var=zeros(nendo,1);
                for i=1:nendo
                    if isempty(find(i_var==i))
                        m_var(i)=i;
                    end
                end
                m_var=nonzeros(m_var);
                lq_instruments.m_var=m_var;
            else
                error('WARNING: There are no instrumnets for ACES!');
            end
        else %if(options_.ACES_solver==1)
            error('WARNING: There are no instrumnets for ACES!');
        end
    end

    % find size xlen of the state vector Y and of A0, A1 and A2 transition matrices:
    % it is the sum the all i variables's lag/lead representations,
    % for each variable i representation being defined as:
    % Max (i_lags-1,0)+ Max (i_leads-1,0)+1
    % so that if variable x appears with 2 lags and 1 lead, and z
    % with 2 lags and 3 leads, the size of the state space is:
    % 1+0+1   +   1+2+1   =6
    % e.g. E_t Y(t+1)=
    %     E_t x(t)
    %     E_t x(t+1)
    %     E_t z(t)
    %     E_t z(t+1)
    %     E_t z(t+2)
    %     E_t z(t+3) 
    
    % partition jacobian:
    jlen=dr.nspred+dr.nsfwrd+M_.endo_nbr+M_.exo_nbr; % length of jacobian
    PSI=-jacobia_(:, jlen-M_.exo_nbr+1:end); % exog
                                             % first transpose M_.lead_lag_incidence';
    lead_lag=M_.lead_lag_incidence';
    max_lead_lag=zeros(nendo,2); % lead/lag representation in Y for each endogenous variable i
    if ( M_.maximum_lag <= 1) && (M_.maximum_lead <= 1)
        xlen=size(jacobia_,1);%nendo;
        AA0=zeros(xlen,xlen);  % empty A0
        AA2=AA0; % empty A2 and A3
        AA3=AA0;
        if xlen==nendo % && M_.maximum_lag <=1 && M_.maximum_lead <=1 % apply a shortcut
            AA1=jacobia_(:,npred+1:npred+nendo); 
            if M_.maximum_lead ==1
                fnd = find(lead_lag(:,M_.maximum_lag+2));
                AA0(:, fnd)= jacobia_(:,nonzeros(lead_lag(:,M_.maximum_lag+2))); %forwd jacobian
            end
            if npred>0 && M_.maximum_lag ==1  
                fnd = find(lead_lag(:,1));
                AA2(:, fnd)= jacobia_(:,nonzeros(lead_lag(:,1))); %backward
            end
        elseif options_.ACES_solver==1 % more endo vars than equations in jacobia_
        if nendo-xlen==num_inst
            PSI=[PSI;zeros(num_inst, M_.exo_nbr)];
            % AA1 contemporary
            AA_all=jacobia_(:,npred+1:npred+nendo); 
            AA1=AA_all(:,lq_instruments.m_var); % endo without instruments
            lq_instruments.ij1=AA_all(:,lq_instruments.inst_var_indices); %  instruments only
            lq_instruments.B1=-[lq_instruments.ij1; eye(num_inst)];
            AA1=[AA1, zeros(xlen,num_inst); zeros(num_inst,xlen), eye(num_inst)];
            %PSI=[PSI; zeros(num_inst,M_.exo_nbr)];
            if M_.maximum_lead ==1 % AA0 forward looking
                AA_all(:,:)=0.0;
                fnd = find(lead_lag(:,M_.maximum_lag+2));
                AA_all(:, fnd)= jacobia_(:,nonzeros(lead_lag(:,M_.maximum_lag+2))); %forwd jacobian
                AA0=AA_all(:,lq_instruments.m_var);
                lq_instruments.ij0=AA_all(:,lq_instruments.inst_var_indices); %  instruments only
                lq_instruments.B0=[lq_instruments.ij0; eye(num_inst)];
                AA0=[AA0, zeros(xlen,num_inst); zeros(num_inst,xlen+num_inst)];
            end
            if npred>0 && M_.maximum_lag ==1  
                AA_all(:,:)=0.0;
                fnd = find(lead_lag(:,1));
                AA_all(:, fnd)= jacobia_(:,nonzeros(lead_lag(:,1))); %backward
                AA2=AA_all(:,lq_instruments.m_var);
                lq_instruments.ij2=AA_all(:,lq_instruments.inst_var_indices); %  instruments only
                lq_instruments.B2=[lq_instruments.ij2; eye(num_inst)];
                AA2=[AA2, lq_instruments.ij2 ; zeros(num_inst,xlen+num_inst)];
            end
        else
            error('ACES number of instruments does match');
        end
        else
            error('More than one lead or lag in the jabian');
        end
        if M_.orig_endo_nbr<nendo
            % findif there are any expecatations at time t
            exp_0= strmatch('AUX_EXPECT_LEAD_0_', M_.endo_names);
            num_exp_0=size(exp_0,1);
            if num_exp_0>0
                AA3(:,exp_0)=AA1(:,exp_0);
                XX0=zeros(nendo,num_exp_0);
                AA1(:,exp_0)=XX0(:,[1:num_exp_0])
            end
        end
    end
    PSI=-[[zeros(xlen-nendo,M_.exo_nbr)];[jacobia_(:, jlen-M_.exo_nbr+1:end)]]; % exog
    cc=0;
    NX=M_.exo_nbr; % no of exogenous varexo shock variables.
    NETA=nfwrd+nboth; % total no of exp. errors  set to no of forward looking equations
    FL_RANK=rank(AA0); % nfwrd+nboth; % min total no of forward looking equations and vars
    
    try
        % call [G1pi,C,impact,nmat,TT1,TT2,gev,eu]=PI_gensys(a0,a1,a2,c,PSI,NX,NETA,NO_FL_EQS)
        % System given as
        %        a0*E_t[y(t+1])+a1*y(t)=a2*y(t-1)+c+psi*eps(t)
        % with eps an exogenous variable process.
        % Returned system is
        %       [s(t)' x(t)' E_t x(t+1)']'=G1pi [s(t-1)' x(t-1)' x(t)]'+C+impact*eps(t),
        %  and (a) the matrix nmat satisfying   nmat*E_t z(t)+ E_t x(t+1)=0
        %      (b) matrices TT1, TT2  that relate y(t) to these states:
        %      y(t)=[TT1 TT2][s(t)' x(t)']'.

        if(options_.ACES_solver==1)
            if isfield(lq_instruments,'xsopt_SS')
                SSbar= diag([lq_instruments.xsopt_SS(m_var)]);% lq_instruments.xsopt_SS(lq_instruments.inst_var_indices)]);
                insSSbar=repmat(lq_instruments.xsopt_SS(lq_instruments.inst_var_indices)',nendo-num_inst,1);
            else
                SSbar= diag([dr.ys(m_var)]);%; dr.ys(lq_instruments.inst_var_indices)]);%(oo_.steady_state);
                insSSbar=repmat(dr.ys(lq_instruments.inst_var_indices)',nendo-num_inst,1);
            end
            SSbar=diag([diag(SSbar);diag(eye(num_inst))]);
            insSSbar=[insSSbar;diag(eye(num_inst))];

            AA0=AA0*SSbar;
            AA1=AA1*SSbar;
            AA2=AA2*SSbar;
            lq_instruments.B1=(lq_instruments.B1).*insSSbar;
        end
        %% for expectational models when complete
        if any(AA3)
            AA3=AA3*SSbar;
            [G1pi,CC,impact,nmat,TT1,TT2,gev,eu, DD, E2,E5, GAMMA, FL_RANK]=PI_gensysEXP(AA0,AA1,-AA2,AA3,cc,PSI,NX,NETA,FL_RANK, M_, options_);
        else
            [G1pi,CC,impact,nmat,TT1,TT2,gev,eu, DD, E2,E5, GAMMA, FL_RANK]=PI_gensys(AA0,AA1,-AA2,AA3,cc,PSI,NX,NETA,FL_RANK, M_, options_);
        end

        % reuse some of the bypassed code and tests that may be needed 
        if (eu(1) ~= 1 || eu(2) ~= 1) && options_.ACES_solver==0
            info(1) = abs(eu(1)+eu(2));
            info(2) = 1.0e+8;
            %     return
        end
        
        dr.PI_ghx=G1pi;
        dr.PI_ghu=impact;
        dr.PI_TT1=TT1;
        dr.PI_TT2=TT2;
        dr.PI_nmat=nmat;
        dr.PI_CC=CC;
        dr.PI_gev=gev;
        dr.PI_eu=eu;
        dr.PI_FL_RANK=FL_RANK;
        %dr.ys=zeros(nendo); % zero steady state
        dr.ghx=G1pi;
        dr.ghu=impact;
        dr.eigval = eig(G1pi);
        dr.rank=FL_RANK;
        
        if options_.ACES_solver==1
            betap=options_.planner_discount;
            sigma_cov=M_.Sigma_e;
            % get W - BY
            W=(1-betap)*GAMMA'*DYN_Q*GAMMA;
            %W=[0]
            ACES.A=G1pi;
            ACES.C=impact; % (:,1);
            ACES.D=DD; %=impact (:,20);
            ACES.E2=E2;
            ACES.E5=E5;
            ACES.GAMMA=GAMMA;
            ACES_M=size(G1pi,2)-FL_RANK;
            ACES_NM=FL_RANK;
            ACES.M=ACES_M;
            ACES.NM=FL_RANK;
            % added by BY
            ACES.Q=DYN_Q;
            ACES.W=W;
            NY=nendo-num_inst;
            
            % save the followings in a subdirectory - BY    
            save ([ACES_DirectoryName,'/',M_.fname '_ACESLQ_Matrices'], 'ACES');
            save ([ACES_DirectoryName,'/',M_.fname '_ACESLQ_GAMMA'], 'GAMMA');
            save ([ACES_DirectoryName,'/',M_.fname '_ACESLQ_A.txt'], 'G1pi', '-ascii', '-double', '-tabs');
            save ([ACES_DirectoryName,'/',M_.fname '_ACESLQ_C.txt'], 'impact','-ascii', '-double', '-tabs');
            save ([ACES_DirectoryName,'/',M_.fname '_ACESLQ_D.txt'], 'DD', '-ascii', '-double', '-tabs');
            save ([ACES_DirectoryName,'/',M_.fname '_ACESLQ_E2.txt'], 'E2', '-ascii', '-double', '-tabs');
            save ([ACES_DirectoryName,'/',M_.fname '_ACESLQ_E5.txt'], 'E5', '-ascii', '-double', '-tabs');
            save ([ACES_DirectoryName,'/',M_.fname '_ACESLQ_GAMMA.txt'], 'GAMMA', '-ascii', '-double', '-tabs');
            %save ([M_.fname '_ACESLQ_M.txt'], 'ACES_M', '-ascii', '-tabs');
            %save ([M_.fname '_ACESLQ_NM.txt'], 'ACES_NM', '-ascii', '-tabs');
            %save ([M_.fname '_ACESLQ_betap.txt'], 'betap', '-ascii', '-tabs');
            %save ([M_.fname '_ACESLQ_NI.txt'], 'num_inst', '-ascii', '-tabs');
            %save ([M_.fname '_ACESLQ_ND.txt'], 'NX', '-ascii', '-tabs');
            %save ([M_.fname '_ACESLQ_NY.txt'], 'NY', '-ascii', '-tabs');
            ACES_VARS=M_.endo_names;
            save ([ACES_DirectoryName,'/',M_.fname '_ACESLQ_VARS.txt'], 'ACES_VARS', '-ascii', '-tabs');
            % added by BY
            % save the char array ACES_VARS into .txt as it is
            fid = fopen(strcat(ACES_DirectoryName,'/',M_.fname,'_ACESLQ_VARnames.txt'),'wt');
            ACES_VARS =[ACES_VARS repmat(sprintf('\n'),size(ACES_VARS,1),1)];
            fwrite(fid,ACES_VARS.');
            fclose(fid);
            % save as integers
            fid = fopen(strcat(ACES_DirectoryName,'/',M_.fname,'_ACESLQ_M.txt'),'wt');
            fprintf(fid,'%d\n',ACES_M);
            fclose(fid);
            fid = fopen(strcat(ACES_DirectoryName,'/',M_.fname,'_ACESLQ_NM.txt'),'wt');
            fprintf(fid,'%d\n',ACES_NM);
            fclose(fid);
            fid = fopen(strcat(ACES_DirectoryName,'/',M_.fname,'_ACESLQ_betap.txt'),'wt');
            fprintf(fid,'%d\n',betap);
            fclose(fid);
            fid = fopen(strcat(ACES_DirectoryName,'/',M_.fname,'_ACESLQ_NI.txt'),'wt');
            fprintf(fid,'%d\n',num_inst);
            fclose(fid);
            fid = fopen(strcat(ACES_DirectoryName,'/',M_.fname,'_ACESLQ_ND.txt'),'wt');
            fprintf(fid,'%d\n',NX);
            fclose(fid);
            fid = fopen(strcat(ACES_DirectoryName,'/',M_.fname,'_ACESLQ_NY.txt'),'wt');
            fprintf(fid,'%d\n',NY);
            fclose(fid);
            save ([ACES_DirectoryName,'/',M_.fname '_ACESLQ_Q.txt'], 'DYN_Q', '-ascii', '-double', '-tabs');
            save ([ACES_DirectoryName,'/',M_.fname '_ACESLQ_W.txt'], 'W', '-ascii', '-double', '-tabs');
            save ([ACES_DirectoryName,'/',M_.fname '_ACESLQ_SIGMAE.txt'], 'sigma_cov', '-ascii', '-double', '-tabs');
        end

    catch
        lerror=lasterror;
        if options_.ACES_solver==1
            disp('Problem with using Part Info ACES solver:');
            error(lerror.message);
        else
            disp('Problem with using Part Info solver');
            error(lerror.message);
        end
    end

    % TODO: 
    % if options_.loglinear == 1
    % if exogenous deterministic variables
    
    return;
