function [dr,info,M_,options_,oo_] = dr_block(dr,task,M_,options_,oo_)
% function [dr,info,M_,options_,oo_] = dr_block(dr,task,M_,options_,oo_)
% computes the reduced form solution of a rational expectation model (first
% approximation of the stochastic model around the deterministic steady state). 
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
%                                 info=6: The jacobian matrix evaluated at the steady state is complex.        
%   M_         [matlab structure]            
%   options_   [matlab structure]
%   oo_        [matlab structure]
%  
% ALGORITHM
%   first order block relaxation method applied to the model
%    E[A Yt-1 + B Yt + C Yt-1 + ut] = 0
%    
% SPECIAL REQUIREMENTS
%   none.
%  

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

info = 0;
verbose = 0;
%it_ = M_.maximum_lag + 1;
z = repmat(dr.ys,1,M_.maximum_lead + M_.maximum_lag + 1);
if (isfield(M_,'block_structure'))
    data = M_.block_structure.block;
    Size = length(M_.block_structure.block);
else
    data = M_;
    Size = 1;
end;
if (options_.bytecode)
    [chck, zz, data]= bytecode('dynamic','evaluate',z,[oo_.exo_simul ...
                        oo_.exo_det_simul], M_.params, dr.ys, 1, data);
else
    [r, data] = feval([M_.fname '_dynamic'], z', [oo_.exo_simul ...
                        oo_.exo_det_simul], M_.params, dr.ys, 2, data);
    chck = 0;
end;
mexErrCheck('bytecode', chck);
dr.rank = 0;
dr.eigval = [];
dr.nstatic = 0;
dr.nfwrd = 0;
dr.npred = 0;
dr.nboth = 0;
dr.nd = 0;
dr.state_var = [];
dr.exo_var = [];
dr.ghx = [];
dr.ghu = [];
for i = 1:Size;
    ghx = [];
    indexi_0 = 0;
    if (verbose)
        disp(['Block ' int2str(i)]);
        disp('-----------');
        data(i)
    end;
    n_pred = data(i).n_backward;
    n_fwrd = data(i).n_forward;
    n_both = data(i).n_mixed;
    n_static = data(i).n_static;
    dr.nstatic = dr.nstatic + n_static;
    dr.nfwrd = dr.nfwrd + n_fwrd;
    dr.npred = dr.npred + n_pred;
    dr.nboth = dr.nboth + n_both;
    nd = n_pred + n_fwrd + 2*n_both;
    dr.nd = dr.nd + nd;
    n_dynamic = n_pred + n_fwrd + n_both;
    exo_nbr = M_.block_structure.block(i).exo_nbr;
    exo_det_nbr = M_.block_structure.block(i).exo_det_nbr;
    jacob = full(data(i).g1);
    lead_lag_incidence = data(i).lead_lag_incidence;
    endo = data(i).variable;
    exo = data(i).exogenous;
    if (verbose)
        disp('jacob');
        disp(jacob);
        disp('lead_lag_incidence');
        disp(lead_lag_incidence);
    end;
    maximum_lag = data(i).maximum_endo_lag;
    maximum_lead = data(i).maximum_endo_lead;
    n = n_dynamic + n_static;
    

    switch M_.block_structure.block(i).Simulation_Type
      case 1
        %Evaluate Forward
        if maximum_lag > 0 && n_pred > 0
            indx_r = find(M_.block_structure.block(i).lead_lag_incidence(1,:));
            indx_c = M_.block_structure.block(i).lead_lag_incidence(1,indx_r);
            data(i).eigval = diag(jacob(indx_r, indx_c));
            data(i).rank = sum(abs(data(i).eigval) > 0);
        else
            data(i).eigval = [];
            data(i).rank = 0;
        end
        dr.eigval = [dr.eigval ; data(i).eigval];
        %First order approximation
        if task ~= 1
            if (maximum_lag > 0)
                indexi_0 = min(lead_lag_incidence(2,:));
                indx_r = find(M_.block_structure.block(i).lead_lag_incidence(1,:));
                indx_c = M_.block_structure.block(i).lead_lag_incidence(1,indx_r);
                ghx = jacob(indx_r, indx_c);
            end;
            ghu = data(i).g1_x;
        end
      case 2
        %Evaluate Backward
        if maximum_lead > 0 && n_fwrd > 0
            indx_r = find(M_.block_structure.block(i).lead_lag_incidence(2,:));
            indx_c = M_.block_structure.block(i).lead_lag_incidence(2,indx_r);
            data(i).eigval = 1./ diag(jacob(indx_r, indx_c));
            data(i).rank = sum(abs(data(i).eigval) > 0);
        else
            data(i).eigval = [];
            data(i).rank = 0;
        end
        dr.rank = dr.rank + data(i).rank;
        dr.eigval = [dr.eigval ; data(i).eigval];
      case 3
        %Solve Forward simple
        if maximum_lag > 0 && n_pred > 0
            data(i).eigval = - jacob(1 , 1 : n_pred) / jacob(1 , n_pred + n_static + 1 : n_pred + n_static + n_pred + n_both);
            data(i).rank = sum(abs(data(i).eigval) > 0);
        else
            data(i).eigval = [];
            data(i).rank = 0;
        end;
        dr.eigval = [dr.eigval ; data(i).eigval];
        %First order approximation
        if task ~= 1
            if (maximum_lag > 0)
                indexi_0 = min(lead_lag_incidence(2,:));
                ghx = - jacob(1 , 1 : n_pred) / jacob(1 , n_pred + n_static + 1 : n_pred + n_static + n_pred + n_both);
            end;
            ghu = data(i).g1_x;
        end
      case 4
        %Solve Backward simple
        if maximum_lead > 0 && n_fwrd > 0
            data(i).eigval = - jacob(1 , n_pred + n - n_fwrd + 1 : n_pred + n) / jacob(1 , n_pred + n + 1 : n_pred + n + n_fwrd) ;
            data(i).rank = sum(abs(data(i).eigval) > 0);
        else
            data(i).eigval = [];
            data(i).rank = 0;
        end;
        dr.rank = dr.rank + data(i).rank;
        dr.eigval = [dr.eigval ; data(i).eigval];
      case 5
        %Solve Forward complete
        if maximum_lag > 0 && n_pred > 0
            data(i).eigval = eig(- jacob(: , 1 : n_pred) / ...
                                 jacob(: , (n_pred + n_static + 1 : n_pred + n_static + n_pred )));
            data(i).rank = sum(abs(data(i).eigval) > 0);
        else
            data(i).eigval = [];
            data(i).rank = 0;
        end;
        dr.eigval = [dr.eigval ; data(i).eigval];
      case 6
        %Solve Backward complete
        if maximum_lead > 0 && n_fwrd > 0
            data(i).eigval = eig(- jacob(: , n_pred + n - n_fwrd + 1: n_pred + n))/ ...
                jacob(: , n_pred + n + 1 : n_pred + n + n_fwrd);
            data(i).rank = sum(abs(data(i).eigval) > 0);
        else
            data(i).eigval = [];
            data(i).rank = 0;
        end;
        dr.rank = dr.rank + data(i).rank;
        dr.eigval = [dr.eigval ; data(i).eigval];
      case 8
        %The lead_lag_incidence contains columns in the following order :
        %  static variables, backward variable, mixed variables and forward variables
        %  
        %Procedes to a QR decomposition on the jacobian matrix to reduce the problem size
        index_c = lead_lag_incidence(2,:);             % Index of all endogenous variables present at time=t
        index_s = lead_lag_incidence(2,1:n_static);    % Index of all static endogenous variables present at time=t
        if n_static > 0
            [Q, R] = qr(jacob(:,index_s));
            aa = Q'*jacob;
        else
            aa = jacob;
        end;
        indexi_0 = min(lead_lag_incidence(2,:));
        index_0m = (n_static+1:n_static+n_pred) + indexi_0 - 1;
        index_0p = (n_static+n_pred+1:n) + indexi_0 - 1;
        index_m = 1:(n_pred+n_both);
        indexi_p = max(lead_lag_incidence(2,:))+1;
        index_p = indexi_p:size(jacob, 2);
        nyf = n_fwrd + n_both;

        A = aa(:,index_m);  % Jacobain matrix for lagged endogeneous variables
        B = aa(:,index_c);  % Jacobian matrix for contemporaneous endogeneous variables
        C = aa(:,index_p);  % Jacobain matrix for led endogeneous variables

        row_indx = n_static+1:n;
        
        D = [[aa(row_indx,index_0m) zeros(n_dynamic,n_both) aa(row_indx,index_p)] ; [zeros(n_both, n_pred) eye(n_both) zeros(n_both, n_both + n_fwrd)]];
        E = [-aa(row_indx,[index_m index_0p])  ; [zeros(n_both, n_both + n_pred) eye(n_both, n_both + n_fwrd) ] ];

        [err, ss, tt, w, sdim, data(i).eigval, info1] = mjdgges(E,D,options_.qz_criterium);
        if (verbose)
            disp('eigval');
            disp(data(i).eigval);
        end;
        if info1
            info(1) = 2;
            info(2) = info1;
            return
        end
        %sdim
        nba = nd-sdim;
        if task == 1
            data(i).rank = rank(w(nd-nyf+1:end,nd-nyf+1:end));
            dr.rank = dr.rank + data(i).rank;
            if ~exist('OCTAVE_VERSION','builtin')
                data(i).eigval = eig(E,D);
            end
            dr.eigval = [dr.eigval ; data(i).eigval];
        end
        if (verbose)
            disp(['sum eigval > 1 = ' int2str(sum(abs(data(i).eigval) > 1.)) ' nyf=' int2str(nyf) ' and dr.rank=' int2str(data(i).rank)]);
            disp(['data(' int2str(i) ').eigval']);
            disp(data(i).eigval);
        end;
        
        %First order approximation
        if task ~= 1
            if nba ~= nyf
                sorted_roots = sort(abs(dr.eigval));
                if isfield(options_,'indeterminacy_continuity')
                    if options_.indeterminacy_msv == 1
                        [ss,tt,w,q] = qz(e',d');
                        [ss,tt,w,q] = reorder(ss,tt,w,q);
                        ss = ss';
                        tt = tt';
                        w  = w';
                        nba = nyf;
                    end
                else
                    if nba > nyf
                        temp = sorted_roots(nd-nba+1:nd-nyf)-1-options_.qz_criterium;
                        info(1) = 3;
                    elseif nba < nyf;
                        temp = sorted_roots(nd-nyf+1:nd-nba)-1-options_.qz_criterium;
                        info(1) = 4;
                    end
                    info(2) = temp'*temp;
                    return
                end
            end
            indx_stable_root = 1: (nd - nyf);    %=> index of stable roots
            indx_explosive_root = n_pred + 1:nd;  %=> index of explosive roots

            % derivatives with respect to dynamic state variables
            % forward variables
            Z = w';
            Z11t = Z(indx_stable_root,    indx_stable_root)';
            Z21  = Z(indx_explosive_root, indx_stable_root);
            Z22  = Z(indx_explosive_root, indx_explosive_root);
            
            if ~isfloat(Z21) && (condest(Z21) > 1e9)
                % condest() fails on a scalar under Octave
                info(1) = 5;
                info(2) = condest(Z21);
                return;
            else
                gx = -inv(Z22) * Z21;
            end
            % predetermined variables
            hx =  Z11t * inv(tt(indx_stable_root, indx_stable_root)) * ss(indx_stable_root, indx_stable_root) * inv(Z11t);
            
            k1 = 1:(n_pred+n_both);
            k2 = 1:(n_fwrd+n_both);

            ghx = [hx(k1,:); gx(k2(n_both+1:end),:)];
            
            if (verbose)
                disp('ghx');
                disp(ghx);
            end;
            
            %lead variables actually present in the model
            
            j4 = n_static+n_pred+1:n_static+n_pred+n_both+n_fwrd;
            j3 = nonzeros(lead_lag_incidence(2,j4)) - n_static - 2 * n_pred - n_both;
            j4 = find(lead_lag_incidence(2,j4));
            if (verbose)
                disp('j3');
                disp(j3);
                disp('j4');
                disp(j4);
            end;

            % derivatives with respect to exogenous variables
            if exo_nbr
                if n_static > 0
                    fu = Q' * data(i).g1_x;
                else
                    fu = data(i).g1_x;
                end;
                
                B_static = [];
                if n_static > 0
                    B_static = B(:,1:n_static);  % submatrix containing the derivatives w.r. to static variables
                end
                B_pred = B(:,n_static+1:n_static+n_pred);
                B_fyd = B(:,n_static+n_pred+1:end);

                ghu = -[B_static C(:,j3)*gx(j4,1:n_pred)+B_pred B_fyd]\fu;
                if (verbose)
                    disp('ghu');
                    disp(ghu);
                end;
            else
                ghu = [];
            end

            % static variables
            if n_static > 0
                temp = - C(1:n_static,j3)*gx(j4,:)*hx;
                if (verbose)
                    disp('temp');
                    disp(temp);
                end;
                j5 = index_m;
                if (verbose)
                    disp('j5');
                    disp(j5);
                end;
                b = aa(:,index_c);
                b10 = b(1:n_static, 1:n_static);
                b11 = b(1:n_static, n_static+1:n);
                if (verbose)
                    disp('b10');
                    disp(b10);
                    disp('b11');
                    disp(b11);
                end;
                temp(:,j5) = temp(:,j5)-A(1:n_static,:);
                if (verbose)
                    disp('temp');
                    disp(temp);
                end;
                disp(temp-b11*ghx);
                temp = b10\(temp-b11*ghx);
                if (verbose)
                    disp('temp');
                    disp(temp);
                end;
                ghx = [temp; ghx];
                temp = [];
                if (verbose)
                    disp('ghx');
                    disp(ghx);
                end;
            end
            
            if options_.loglinear == 1
                k = find(dr.kstate(:,2) <= M_.maximum_endo_lag+1);
                klag = dr.kstate(k,[1 2]);
                k1 = dr.order_var;
                
                ghx = repmat(1./dr.ys(k1),1,size(ghx,2)).*ghx.* ...
                      repmat(dr.ys(k1(klag(:,1)))',size(ghx,1),1);
                ghu = repmat(1./dr.ys(k1),1,size(ghu,2)).*ghu;
            end

            if options_.aim_solver ~= 1 && options_.use_qzdiv
                % Necessary when using Sims' routines for QZ
                gx = real(gx);
                hx = real(hx);
                ghx = real(ghx);
                ghu = real(ghu);
            end
            ghx
            %exogenous deterministic variables
            if exo_det_nbr > 0
                f1 = sparse(jacobia_(:,nonzeros(M_.lead_lag_incidence(M_.maximum_endo_lag+2:end,order_var))));
                f0 = sparse(jacobia_(:,nonzeros(M_.lead_lag_incidence(M_.maximum_endo_lag+1,order_var))));
                fudet = data(i).g1_xd;
                M1 = inv(f0+[zeros(n,n_static) f1*gx zeros(n,nyf-n_both)]);
                M2 = M1*f1;
                dr.ghud = cell(M_.exo_det_length,1);
                dr.ghud{1} = -M1*fudet;
                for i = 2:M_.exo_det_length
                    dr.ghud{i} = -M2*dr.ghud{i-1}(end-nyf+1:end,:);
                end
            end
            
            
            %Endogeneous variables of the previous blocks
            
            
        end
    end;
    if task ~=1
        if (maximum_lag > 0)
            lead_lag_incidence(maximum_lag+1, n_static+1:n_static + n_pred + n_both) - indexi_0 + 1
            state_var = endo(lead_lag_incidence(maximum_lag+1, n_static+1:n_static + n_pred + n_both) - indexi_0 + 1);
            [common_state_var, indx_common_dr_state_var, indx_common_state_var]  = intersect(dr.state_var, state_var);
            [diff_state_var, indx_diff_dr_state_var, indx_diff_state_var] =  setxor(dr.state_var, state_var);
            [union_state_var, indx_union_dr_state_var, indx_union_state_var] = union(dr.state_var, state_var);
            [row_dr_ghx, col_dr_ghx] = size(dr.ghx);
            ghx_new = zeros(row_dr_ghx + n, length(union_state_var));
            ghx_new(1:row_dr_ghx,  1:col_dr_ghx) = dr.ghx;
            ghx_new(row_dr_ghx + 1: row_dr_ghx + n, indx_common_dr_state_var) = ghx(:, indx_common_state_var);
            ghx_new(row_dr_ghx + 1: row_dr_ghx + n, length(dr.state_var)+1:length(dr.state_var)+length(indx_diff_state_var)) = ghx(:, indx_diff_state_var);
            dr.ghx = ghx_new;
            dr.state_var = [dr.state_var state_var(indx_diff_state_var)];
        end;
        exo_var = exo;
        [common_exo_var, indx_common_dr_exo_var, indx_common_exo_var]  = intersect(dr.exo_var, exo_var);
        [diff_exo_var, indx_diff_dr_exo_var, indx_diff_exo_var] =  setxor(dr.exo_var, exo_var);
        [union_exo_var, indx_union_dr_exo_var, indx_union_exo_var] = union(dr.exo_var, exo_var);
        [row_dr_ghu, col_dr_ghu] = size(dr.ghu);
        ghu_new = zeros(row_dr_ghu + exo_nbr, length(union_exo_var));
        ghu_new(1:row_dr_ghu,  1:col_dr_ghu) = dr.ghu;
        ghu_new(row_dr_ghu + 1: row_dr_ghu + n, indx_common_dr_exo_var) = ghu(:, indx_common_exo_var);
        ghu_new(row_dr_ghu + 1: row_dr_ghu + n, length(dr.exo_var)+1:length(dr.exo_var)+length(indx_diff_exo_var)) = ghu(:, indx_diff_exo_var);
        dr.ghu = ghu_new;
        dr.exo_var = [dr.exo_var exo_var(indx_diff_exo_var)];
    end
end;
if (verbose)
    dr.ghx
    dr.ghu
end;
if (task == 1)
    return;
end;