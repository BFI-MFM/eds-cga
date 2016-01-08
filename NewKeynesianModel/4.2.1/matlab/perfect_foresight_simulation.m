function info = perfect_foresight_simulation(compute_linear_solution,steady_state)
% Performs deterministic simulations with lead or lag on one period
%
% INPUTS
%   endo_simul                  [double]     n*T matrix, where n is the number of endogenous variables.
%   exo_simul                   [double]     q*T matrix, where q is the number of shocks.
%   compute_linear_solution     [integer]    scalar equal to zero or one.     
%     
% OUTPUTS
%   none
%
% ALGORITHM
%   Laffargue, Boucekkine, Juillard (LBJ)
%   see Juillard (1996) Dynare: A program for the resolution and
%   simulation of dynamic models with forward variables through the use
%   of a relaxation algorithm. CEPREMAP. Couverture Orange. 9602.
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 2009-2010 Dynare Team
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

global M_ options_ it_ oo_

persistent lead_lag_incidence dynamic_model ny nyp nyf nrs nrc iyf iyp isp is isf isf1 iz icf ghx iflag

if ~nargin && isempty(iflag)% Initialization of the persistent variables.
    lead_lag_incidence = M_.lead_lag_incidence; 
    dynamic_model = [M_.fname '_dynamic'];
    ny   = size(oo_.endo_simul,1); 
    nyp  = nnz(lead_lag_incidence(1,:));% number of lagged variables.
    nyf  = nnz(lead_lag_incidence(3,:));% number of leaded variables. 
    nrs  = ny+nyp+nyf+1;
    nrc  = nyf+1; 
    iyf  = find(lead_lag_incidence(3,:)>0);% indices for leaded variables. 
    iyp  = find(lead_lag_incidence(1,:)>0);% indices for lagged variables. 
    isp  = 1:nyp;
    is   = (nyp+1):(nyp+ny); % Indices for contemporaneaous variables. 
    isf  = iyf+nyp; 
    isf1 = (nyp+ny+1):(nyf+nyp+ny+1);     
    iz   = 1:(ny+nyp+nyf);
    icf  = 1:size(iyf,2);
    info = [];
    iflag = 1;
    return
end

initial_endo_simul = oo_.endo_simul;

warning_old_state = warning;
warning off all

if nargin<1
    compute_linear_solution = 0;
else
    if nargin<2
        error('The steady state (second input argument) is missing!');
    end
end

if ~isstruct(compute_linear_solution) && compute_linear_solution 
    [dr,info]=resol(steady_state,0);
elseif isstruct(compute_linear_solution)
    dr = compute_linear_solution;
    compute_linear_solution = 1;
end

if compute_linear_solution
    ghx(dr.order_var,:) = dr.ghx;
    ghx = ghx(iyf,:);
end

periods = options_.periods; 

stop    = 0 ; 
it_init = M_.maximum_lag+1; 

info.convergence = 1; 
info.time  = 0; 
info.error = 0; 
info.iterations.time  = zeros(options_.maxit_,1); 
info.iterations.error = info.iterations.time; 

last_line = options_.maxit_;
error_growth = 0;

h1 = clock;
for iter = 1:options_.maxit_ 
    h2 = clock;
    if options_.terminal_condition
        c = zeros(ny*(periods+1),nrc);
    else
        c = zeros(ny*periods,nrc);
    end
    it_ = it_init;
    z = [ oo_.endo_simul(iyp,it_-1) ; oo_.endo_simul(:,it_) ; oo_.endo_simul(iyf,it_+1) ]; 
    [d1,jacobian] = feval(dynamic_model,z,oo_.exo_simul, M_.params, it_); 
    jacobian = [jacobian(:,iz) , -d1]; 
    ic = 1:ny;
    icp = iyp;
    c(ic,:) = jacobian(:,is)\jacobian(:,isf1) ; 
    for it_ = it_init+(1:periods-1-(options_.terminal_condition==2))
        z = [ oo_.endo_simul(iyp,it_-1) ; oo_.endo_simul(:,it_) ; oo_.endo_simul(iyf,it_+1)]; 
        [d1,jacobian] = feval(dynamic_model,z,oo_.exo_simul, M_.params, it_); 
        jacobian = [jacobian(:,iz) , -d1];
        jacobian(:,[isf nrs]) = jacobian(:,[isf nrs])-jacobian(:,isp)*c(icp,:); 
        ic = ic + ny;
        icp = icp + ny;
        c(ic,:) = jacobian(:,is)\jacobian(:,isf1); 
    end
    if options_.terminal_condition
        if options_.terminal_condition==1% Terminal condition is Y_{T} = Y_{T+1} 
            s = eye(ny);
            s(:,isf) = s(:,isf)+c(ic,1:nyf);
            ic = ic + ny;
            c(ic,nrc) = s\c(ic,nrc);
        else% Terminal condition is Y_{T+1}-Y^{\star} = TransitionMatrix*(Y_{T}-Y^{\star})
            it_ = it_+1;
            z = [ oo_.endo_simul(iyp,it_-1) ; oo_.endo_simul(:,it_) ; oo_.endo_simul(iyf,it_+1) ] ;
            [d1,jacobian] = feval(dynamic_model,z,oo_.exo_simul, M_.params, it_);
            jacobian = [jacobian(:,iz) -d1];
            jacobian(:,[isf nrs]) = jacobian(:,[isf nrs])-jacobian(:,isp)*c(icp,:) ;
            ic = ic + ny;
            icp = icp + ny;
            s = jacobian(:,is);
            s(:,iyp) = s(:,iyp)+jacobian(:,isf)*ghx;
            c (ic,:) = s\jacobian(:,isf1);
        end
        c = bksup0(c,ny,nrc,iyf,icf,periods);
        c = reshape(c,ny,periods+1);
        oo_.endo_simul(:,it_init+(0:periods)) = oo_.endo_simul(:,it_init+(0:periods))+options_.slowc*c;
    else% Terminal condition is Y_{T}=Y^{\star}
        c = bksup0(c,ny,nrc,iyf,icf,periods);
        c = reshape(c,ny,periods);
        oo_.endo_simul(:,it_init+(0:periods-1)) = oo_.endo_simul(:,it_init+(0:periods-1))+options_.slowc*c; 
    end
    err = max(max(abs(c))); 
    info.iterations.time(iter)  = etime(clock,h2); 
    info.iterations.error(iter) = err;
    if iter>1
        error_growth = error_growth + (info.iterations.error(iter)>info.iterations.error(iter-1));
    end
    if isnan(err) || error_growth>3
        last_line = iter;
        break
    end
    if err < options_.dynatol
        stop = 1;
        info.time  = etime(clock,h1); 
        info.error = err;
        info.iterations.time = info.iterations.time(1:iter); 
        info.iterations.error  = info.iterations.error(1:iter);
        break
    end
end

if stop && options_.terminal_condition==2
    % Compute the distance to the deterministic steady state (for the subset of endogenous variables with a non zero 
    % steady state) at the last perdiod.
    idx = find(abs(oo_.steady_state)>0);
    distance_to_steady_state = abs(((oo_.endo_simul(idx,end)-oo_.steady_state(idx))./oo_.steady_state(idx)))*100;
    disp(['(max) Distance to steady state at the end is (in percentage):' num2str(max(distance_to_steady_state))])
end

if ~stop
    info.time  = etime(clock,h1);
    info.error = err;
    info.convergence = 0;
    info.iterations.time  = info.iterations.time(1:last_line);
    info.iterations.error = info.iterations.error(1:last_line);
    oo_.endo_simul = initial_endo_simul;
end

warning(warning_old_state);
