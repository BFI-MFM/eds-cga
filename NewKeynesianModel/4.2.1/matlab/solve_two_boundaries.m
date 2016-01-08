function y = solve_two_boundaries(fname, y, x, params, steady_state, y_index, nze, periods, y_kmin_l, y_kmax_l, is_linear, Block_Num, y_kmin, maxit_, solve_tolf, lambda, cutoff, stack_solve_algo)
% Computes the deterministic simulation of a block of equation containing
% both lead and lag variables using relaxation methods 
%
% INPUTS
%   fname               [string]        name of the file containing the block
%                                       to simulate
%   y                   [matrix]        All the endogenous variables of the model
%   x                   [matrix]        All the exogenous variables of the model
%   params              [vector]        All the parameters of the model
%   steady_state        [vector]        steady state of the model
%   y_index             [vector of int] The index of the endogenous variables of
%                                       the block
%   nze                 [integer]       number of non-zero elements in the
%                                       jacobian matrix
%   periods             [integer]       number of simulation periods
%   y_kmin_l            [integer]       maximum number of lag in the block
%   y_kmax_l            [integer]       maximum number of lead in the block
%   is_linear           [integer]       if is_linear=1 the block is linear
%                                       if is_linear=0 the block is not linear
%   Block_Num           [integer]       block number
%   y_kmin              [integer]       maximum number of lag in the model
%   maxit_              [integer]       maximum number of iteration in Newton
%   solve_tolf          [double]        convergence criteria
%   lambda              [double]        initial value of step size in
%   Newton
%   cutoff              [double]        cutoff to correct the direction in Newton in case
%                                       of singular jacobian matrix
%   stack_solve_algo    [integer]       linear solver method used in the
%                                       Newton algorithm : 
%                                            - 1 sprse LU
%                                            - 2 GMRES
%                                            - 3 BicGStab
%                                            - 4 Optimal path length
%
% OUTPUTS
%   y                   [matrix]        All endogenous variables of the model      
%  
% ALGORITHM
%   Newton with LU or GMRES or BicGstab
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

global oo_ M_;
cvg=0;
iter=0;
Per_u_=0;
g2 = [];
g3 = [];
Blck_size=size(y_index,2);
correcting_factor=0.01;
luinc_tol=1e-10;
max_resa=1e100;
Jacobian_Size=Blck_size*(y_kmin+y_kmax_l +periods);
g1=spalloc( Blck_size*periods, Jacobian_Size, nze*periods);
reduced = 0;
while ~(cvg==1 || iter>maxit_),
    [r, y, g1, g2, g3, b]=feval(fname, y, x, params, steady_state, periods, 0, y_kmin, Blck_size);
    %     fjac = zeros(Blck_size, Blck_size*(y_kmin_l+1+y_kmax_l));
    %     disp(['Blck_size=' int2str(Blck_size) ' size(y_index)=' int2str(size(y_index,2))]);
    %     dh = max(abs(y(y_kmin+1-y_kmin_l:y_kmin+1+y_kmax_l, y_index)),options_.gstep*ones(y_kmin_l+1+y_kmax_l, Blck_size))*eps^(1/3);
    %     fvec = r;
    %     %for i = y_kmin+1-y_kmin_l:y_kmin+1+y_kmax_l
    %     i = y_kmin+1;
    %       i
    %       for j = 1:Blck_size
    %             ydh = y ;
    %           ydh(i, y_index(j)) = ydh(i, y_index(j)) + dh(i, j)  ;
    %           if(j==11 && i==2)
    %               disp(['y(i,y_index(11)=' int2str(y_index(11)) ')= ' num2str(y(i,y_index(11))) ' ydh(i, y_index(j))=' num2str(ydh(i, y_index(j))) ' dh(i,j)= ' num2str(dh(i,j))]);
    %           end;
    %           [t, y1, g11, g21, g31, b1]=feval(fname, ydh, x, params, periods, 0, y_kmin, Blck_size);
    %           fjac(:,j+(i-(y_kmin+1-y_kmin_l))*Blck_size) = (t(:, 1+y_kmin) - fvec(:, 1+y_kmin))./dh(i, j) ;
    %           if(j==11 && i==2)
    %                disp(['fjac(:,' int2str(j+(i-(y_kmin+1-y_kmin_l))*Blck_size) ')=']);
    %                disp([num2str(fjac(:,j+(i-(y_kmin+1-y_kmin_l))*Blck_size))]);
    %           end;
    %       end;
    % %    end
    %     %diff = g1(1:Blck_size, 1:Blck_size*(y_kmin_l+1+y_kmax_l)) -fjac;
    %     diff = g1(1:Blck_size, y_kmin_l*Blck_size+1:(y_kmin_l+1)*Blck_size) -fjac(1:Blck_size, y_kmin_l*Blck_size+1:(y_kmin_l+1)*Blck_size);
    %     disp(diff);
    %     [c_max, i_c_max] = max(abs(diff));
    %     [l_c_max, i_r_max] = max(c_max);
    %     disp(['maximum element row=' int2str(i_c_max(i_r_max)) ' and column=' int2str(i_r_max) ' value = ' num2str(l_c_max)]);
    %     equation = i_c_max(i_r_max);
    %     variable = i_r_max;
    %     variable
    %     disp(['equation ' int2str(equation) ' and variable ' int2str(y_index(mod(variable, Blck_size))) ' ' M_.endo_names(y_index(mod(variable, Blck_size)), :)]);
    %     disp(['g1(' int2str(equation) ', ' int2str(variable) ')=' num2str(g1(equation, y_kmin_l*Blck_size+variable),'%3.10f') ' fjac(' int2str(equation) ', ' int2str(variable) ')=' num2str(fjac(equation, y_kmin_l*Blck_size+variable), '%3.10f')]);
    %     return;



    %     for i=1:periods;
    %       disp([sprintf('%5.14f ',[T9025(i) T1149(i) T11905(i)])]);
    %     end;
    %     return;
    %residual = r(:,y_kmin+1:y_kmin+1+y_kmax_l);
    %num2str(residual,' %1.6f')
    %jac_ = g1(1:(y_kmin)*Blck_size, 1:(y_kmin+1+y_kmax_l)*Blck_size);
    %jac_
    
    g1a=g1(:, y_kmin*Blck_size+1:(periods+y_kmin)*Blck_size);
    term1 = g1(:, 1:y_kmin_l*Blck_size)*reshape(y(1+y_kmin-y_kmin_l:y_kmin,y_index)',1,y_kmin_l*Blck_size)';
    term2 = g1(:, (periods+y_kmin_l)*Blck_size+1:(periods+y_kmin_l+y_kmax_l)*Blck_size)*reshape(y(periods+y_kmin+1:periods+y_kmin+y_kmax_l,y_index)',1,y_kmax_l*Blck_size)';
    b = b - term1 - term2;
    
    %      fid = fopen(['result' num2str(iter)],'w');
    %      fg1a = full(g1a);
    %      fprintf(fid,'%d\n',size(fg1a,1));
    %      fprintf(fid,'%d\n',size(fg1a,2));
    %      fprintf(fid,'%5.14f\n',fg1a);
    %      fprintf(fid,'%d\n',size(b,1));
    %      fprintf(fid,'%5.14f\n',b);
    %      fclose(fid);
    %      return;
    %ipconfigb_ = b(1:(1+y_kmin)*Blck_size);
    %b_ 
    
    
    [max_res, max_indx]=max(max(abs(r')));
    if(~isreal(r))
        max_res = (-max_res^2)^0.5;
    end;
    %     if(~isreal(r))
    %       max_res=(-(max(max(abs(r))))^2)^0.5;
    %     else
    %       max_res=max(max(abs(r)));
    %     end;
    if(~isreal(max_res) || isnan(max_res))
        cvg = 0;
    elseif(is_linear && iter>0)
        cvg = 1;
    else
        cvg=(max_res<solve_tolf);
    end;
    if(~cvg)
        if(iter>0)
            if(~isreal(max_res) || isnan(max_res) || (max_resa<max_res && iter>1))
                if(~isreal(max_res))
                    disp(['Variable ' M_.endo_names(max_indx,:) ' (' int2str(max_indx) ') returns an undefined value']);
                end;
                if(isnan(max_res))
                    detJ=det(g1aa);
                    if(abs(detJ)<1e-7)
                        max_factor=max(max(abs(g1aa)));
                        ze_elem=sum(diag(g1aa)<cutoff);
                        disp([num2str(full(ze_elem),'%d') ' elements on the Jacobian diagonal are below the cutoff (' num2str(cutoff,'%f') ')']);
                        if(correcting_factor<max_factor)
                            correcting_factor=correcting_factor*4;
                            disp(['The Jacobain matrix is singular, det(Jacobian)=' num2str(detJ,'%f') '.']);
                            disp(['    trying to correct the Jacobian matrix:']);
                            disp(['    correcting_factor=' num2str(correcting_factor,'%f') ' max(Jacobian)=' num2str(full(max_factor),'%f')]);
                            dx = (g1aa+correcting_factor*speye(periods*Blck_size))\ba- ya;
                            y(1+y_kmin:periods+y_kmin,y_index)=reshape((ya_save+lambda*dx)',length(y_index),periods)';
                            continue;
                        else
                            disp('The singularity of the jacobian matrix could not be corrected');
                            return;
                        end;
                    end;
                elseif(lambda>1e-8)
                    lambda=lambda/2;
                    reduced = 1;
                    disp(['reducing the path length: lambda=' num2str(lambda,'%f')]);
                    y(1+y_kmin:periods+y_kmin,y_index)=reshape((ya_save+lambda*dx)',length(y_index),periods)';
                    continue;
                else
                    if(cutoff == 0)
                        fprintf('Error in simul: Convergence not achieved in block %d, after %d iterations.\n Increase "options_.maxit_".\n',Block_Num, iter);
                    else
                        fprintf('Error in simul: Convergence not achieved in block %d, after %d iterations.\n Increase "options_.maxit_" or set "cutoff=0" in model options.\n',Block_Num, iter);
                    end;
                    oo_.deterministic_simulation.status = 0;
                    oo_.deterministic_simulation.error = max_res;
                    oo_.deterministic_simulation.iterations = iter;
                    oo_.deterministic_simulation.block(Block_Num).status = 0;% Convergency failed.
                    oo_.deterministic_simulation.block(Block_Num).error = max_res;
                    oo_.deterministic_simulation.block(Block_Num).iterations = iter;
                    return;
                end;
            else
                if(lambda<1)
                    lambda=max(lambda*2, 1);
                end;
            end;
        end;
        ya = reshape(y(y_kmin+1:y_kmin+periods,y_index)',1,periods*Blck_size)';
        ya_save=ya;
        g1aa=g1a;
        ba=b;
        max_resa=max_res;
        if(stack_solve_algo==0),
            dx = g1a\b- ya;
            ya = ya + lambda*dx;
            y(1+y_kmin:periods+y_kmin,y_index)=reshape(ya',length(y_index),periods)';
        elseif(stack_solve_algo==1),
            for t=1:periods;
                first_elem = (t-1)*Blck_size+1;
                last_elem = t*Blck_size;
                next_elem = (t+1)*Blck_size;
                Elem = first_elem:last_elem;
                Elem_1 = last_elem+1:next_elem;
                B1_inv = inv(g1a(Elem, Elem));
                if (t < periods)
                    S1 = B1_inv * g1a(Elem, Elem_1);
                end;
                g1a(Elem, Elem_1) = S1;
                b(Elem) = B1_inv * b(Elem);
                g1a(Elem, Elem) = ones(Blck_size, Blck_size);
                if (t < periods)
                    g1a(Elem_1, Elem_1) = g1a(Elem_1, Elem_1) - g1a(Elem_1, Elem) * S1;
                    b(Elem_1) = b(Elem_1) - g1a(Elem_1, Elem) * b(Elem);
                    g1a(Elem_1, Elem) = zeros(Blck_size, Blck_size);
                end;
            end;
            za = b(Elem);
            zaa = za;
            y_Elem = (periods - 1) * Blck_size + 1:(periods) * Blck_size;
            dx = ya;
            dx(y_Elem) = za - ya(y_Elem);
            ya(y_Elem) = ya(y_Elem) + lambda*dx(y_Elem);
            for t=periods-1:-1:1;
                first_elem = (t-1)*Blck_size+1;
                last_elem = t*Blck_size;
                next_elem = (t+1)*Blck_size;
                Elem_1 = last_elem+1:next_elem;
                Elem = first_elem:last_elem;
                za = b(Elem) - g1a(Elem, Elem_1) * zaa;
                zaa = za;
                y_Elem = Blck_size * (t-1)+1:Blck_size * (t);
                dx(y_Elem) = za - ya(y_Elem);
                ya(y_Elem) = ya(y_Elem) + lambda*dx(y_Elem);
                y(y_kmin + t, y_index) = ya(y_Elem);
            end;
        elseif(stack_solve_algo==2),
            flag1=1;
            while(flag1>0)
                [L1, U1]=luinc(g1a,luinc_tol);
                [za,flag1] = gmres(g1a,b,Blck_size,1e-6,Blck_size*periods,L1,U1);
                if (flag1>0 || reduced)
                    if(flag1==1)
                        disp(['Error in simul: No convergence inside GMRES after ' num2str(periods*10,'%6d') ' iterations, in block' num2str(Block_Size,'%3d')]);
                    elseif(flag1==2)
                        disp(['Error in simul: Preconditioner is ill-conditioned, in block' num2str(Block_Size,'%3d')]);
                    elseif(flag1==3)
                        disp(['Error in simul: GMRES stagnated (Two consecutive iterates were the same.), in block' num2str(Block_Size,'%3d')]);
                    end;
                    luinc_tol = luinc_tol/10;
                    reduced = 0;
                else
                    dx = za - ya;
                    ya = ya + lambda*dx;
                    y(1+y_kmin:periods+y_kmin,y_index)=reshape(ya',length(y_index),periods)';
                end;
            end;
        elseif(stack_solve_algo==3),
            flag1=1;
            while(flag1>0)
                [L1, U1]=luinc(g1a,luinc_tol);
                [za,flag1] = bicgstab(g1a,b,1e-7,Blck_size*periods,L1,U1);
                if (flag1>0 || reduced)
                    if(flag1==1)
                        disp(['Error in simul: No convergence inside BICGSTAB after ' num2str(periods*10,'%6d') ' iterations, in block' num2str(Block_Size,'%3d')]);
                    elseif(flag1==2)
                        disp(['Error in simul: Preconditioner is ill-conditioned, in block' num2str(Block_Size,'%3d')]);
                    elseif(flag1==3)
                        disp(['Error in simul: GMRES stagnated (Two consecutive iterates were the same.), in block' num2str(Block_Size,'%3d')]);
                    end;
                    luinc_tol = luinc_tol/10;
                    reduced = 0;
                else
                    dx = za - ya;
                    ya = ya + lambda*dx;
                    y(1+y_kmin:periods+y_kmin,y_index)=reshape(ya',length(y_index),periods)';
                end;
            end;
        elseif(stack_solve_algo==4),
            ra = reshape(r(:, y_kmin+1:periods+y_kmin),periods*Blck_size, 1);
            stpmx = 100 ;
            stpmax = stpmx*max([sqrt(ya'*ya);size(y_index,2)]);
            nn=1:size(ra,1);
            g = (ra'*g1a)';
            f = 0.5*ra'*ra;
            p = -g1a\ra;
            [yn,f,ra,check]=lnsrch1(ya,f,g,p,stpmax, ...
                                    'lnsrch1_wrapper_two_boundaries',nn,nn,  fname, y, y_index, x, params, steady_state, periods, y_kmin, Blck_size);
            dx = ya - yn;
            y(1+y_kmin:periods+y_kmin,y_index)=reshape(yn',length(y_index),periods)';
        end
    end
    iter=iter+1;
    disp(['iteration: ' num2str(iter,'%d') ' error: ' num2str(max_res,'%e')]);
end;
if (iter>maxit_)
    disp(['No convergence after ' num2str(iter,'%4d') ' iterations in Block ' num2str(Block_Num,'%d')]);
    oo_.deterministic_simulation.status = 0;
    oo_.deterministic_simulation.error = max_res;
    oo_.deterministic_simulation.iterations = iter;
    oo_.deterministic_simulation.block(Block_Num).status = 0;% Convergency failed.
    oo_.deterministic_simulation.block(Block_Num).error = max_res;
    oo_.deterministic_simulation.block(Block_Num).iterations = iter;
    return;
end
oo_.deterministic_simulation.status = 1;
oo_.deterministic_simulation.error = max_res;
oo_.deterministic_simulation.iterations = iter;
oo_.deterministic_simulation.block(Block_Num).status = 1;% Convergency obtained.
oo_.deterministic_simulation.block(Block_Num).error = max_res;
oo_.deterministic_simulation.block(Block_Num).iterations = iter;
return;
