% NK_Degree3.m is a routine for computing degree-three polynomial solutions  
% to a new Keynesian model using an epsilon distinguishable set algorithm 
% (EDS) and cluster grid algorithm (CGA), as described in the article  
% "Merging simulation and projection approaches to solve high-dimensional  
% problems with an application to a new Keynesian model" by Lilia Maliar  
% and Serguei Maliar (2015), Quantitative Economics 6/1, pages 1–47 
% (henceforth, MM, 2015). This routine must be ran after "Main_NK_Degree2.m".
%
%
% This version: March 19, 2015. First version: June 27, 2011.
%
% -------------------------------------------------------------------------
% Copyright © 2015 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% The code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

clc;
clear all;

tic
load pie1_D2;    % Load the second-degree approximation, computed by 
                 % "Main_NK_Degree2.m"

% 1. The number of state variables, simulation length, and polynomial 
% degree
% ---------------------------------------------------------------------
N     = 8;       % Indicate the number of state variables
T     = 10000;   % Choose the simulation length for the solution procedure,
                 % T<=10,000   
Degree = 3;      % Consider the third-degree polynomial approximation


% 2. An initial condition for the polynomial coefficients
%--------------------------------------------------------
X0 = Ord_Polynomial_N(Data(1:T-1,:),Degree);
% Construct the matrix of explanatory variables X0 on the series of state 
% variables from the previously computed time-series perturbation solution; 
% columns of X0 are given by the basis functions of the polynomial of degree 
% "Degree" 

npol_2d = size(X0,2);       % Number of coefficients in polynomial of degree 
                            % "Degree" (recall that we chose Degree=3)
                            
VK = zeros(npol_2d,3,2);    % Initialize the matrix of polynomial coefficients; 
                            % npol_2d-by-3-by-2, where 3 is the number of 
                            % decision rules we compute, i.e., for S, F,  
                            % (C^(-gam) (i.e., marginal utility), and 2 
                            % stands for 2 degrees of approximations we  
                            % consider, of degree 1 and degree 2; npol_2d-by-3-by-2 

VK(1:45,1:3,1)   = vk_2d;   % Fill in VK with the known second-degree 
                            % polynomial coefficients

% 3. Initialize the third-degree polynomial approximations of the policy 
% functions
% ----------------------------------------------------------------------
vk_2d = VK(:,:,1);         % Start from the second-degree polynomial solution
                           % obtained under EDS using "Main_NK_Degree2.m"

                           
tic;                       % Start counting time for computing the EDS solution

% 4. Select the integration method 
% ---------------------------------
[n_nodes,epsi_nodes,weight_nodes] = Monomials_1(6,vcv);
                             % Monomial integration rule with 2N nodes
                             
%[n_nodes,epsi_nodes,weight_nodes] = Monomials_2(N,vcv);
                             % Monomial integration rule with 2N^2+1 nodes

%[n_nodes,epsi_nodes,weight_nodes] = Monomials_3(N,vcv);
                             % Monomial integration rule with 2^N nodes;  
                             % this rule coincides with Gauss-Hermite 
                             % quadrature (product) integration rule with 
                             % 2 nodes in each dimension 

%n_nodes = 1; epsi_nodes = zeros(1,N); weight_nodes = 1;                             
                             % Gauss-Hermite quadrature integration rule 
                             % with one node 


% 5. Compute the number of grid points
%------------------------------------
 n_G = size(Grid,1);    % The number of points in the EDS grid; the obtained 
                        % grid is of size n_G-by-N

% 6. Allocate memory to the integrals in the right side of the 3 Euler 
% equations that we parameterize; see conditions (32), (33), (34) in MM (2015)
%-------------------------------------------------------------------------
 e = zeros(n_G,3);     

% 7. Form complete polynomial of a degree "Degree" on the grid points for  
% future regressions
%-------------------------------------------------------------------------
 X0_G = Ord_Polynomial_N(Grid,Degree);                              

% 8. Allocate memory to the next-period values of S, F and C on the grid 
%------------------------------------------------------------------------
   % Allocate memory to S, F, C obtained in the previous iteration (to check
   % convergence)
   S0_old_G = ones(n_G,1);
   F0_old_G = ones(n_G,1);
   C0_old_G = ones(n_G,1);
   
   % Allocate memory to S, F, C obtained in the current iteration (to check
   % convergence) 
   S0_new_G = ones(n_G,1);     
   F0_new_G = ones(n_G,1);     
   C0_new_G = ones(n_G,1);     
                             
% 9. The parameters of the EDS algorithm 
% -----------------------------------------
 damp     = [0.1 0.1 0.1];% Damping parameter for (fixed-point) iteration on 
                          % the coefficients of the 3 policy functions (for
                          % S, F and C^(-gam))
 dif_EDS	  = 1e+10;    % Convergence criterion (initially is not satisfied)

% 10. The main iterative cycle of the EDS algorithm
% --------------------------------------------------              
 while dif_EDS > 1e-7;% The convergence criterion (which is unit free 
                      % because dif_EDS is unit free)
    
    % 10.1 Computations for each grid point
    %--------------------------------------
    for i = 1:n_G;    % For each grid point, ... 

        % 10.1.1 Variables in a grid point i 
        % ----------------------------------   
        
         % Endogenous state variables on the grid
         %---------------------------------------
         R0  = exp(Grid(i,1));
         delta0  = exp(Grid(i,2));
        
         % Exogenous state variables on the grid
         %--------------------------------------
         nuR0 = Grid(i,3);          % The variable "nuR" appears in the 3d
                                    % column of "Grid"
         nua0 = Grid(i,4);            
         nuL0 = Grid(i,5);           
         nuu0 = Grid(i,6);            
         nuB0 = Grid(i,7);           
         nuG0 = Grid(i,8);            
        
         % Compute future shocks in all grid points and all integration nodes
         % ------------------------------------------------------------------                           
         nuR1(1:n_nodes,1) = (ones(n_nodes,1)*nuR0)*rho_nuR+epsi_nodes(:,1); 
         nua1(1:n_nodes,1) = (ones(n_nodes,1)*nua0)*rho_nua+epsi_nodes(:,2); 
         nuL1(1:n_nodes,1) = (ones(n_nodes,1)*nuL0)*rho_nuL+epsi_nodes(:,3); 
         nuu1(1:n_nodes,1) = (ones(n_nodes,1)*nuu0)*rho_nuu+epsi_nodes(:,4); 
         nuB1(1:n_nodes,1) = (ones(n_nodes,1)*nuB0)*rho_nuB+epsi_nodes(:,5); 
         nuG1(1:n_nodes,1) = (ones(n_nodes,1)*nuG0)*rho_nuG+epsi_nodes(:,6);
         % The size of each od these vectors is 1-by-n_nodes


        % Give a name to an ith row of previously created complete polynomial 
        % "X0_G" (the one that corresponds to the grid point i)
        %------------------------------------------------------------------
        X0 = X0_G(i,:);                 % This is done for convenience 
                                    
        % 10.1.2 Current-period choices in a grid point i
        % -----------------------------------------------
        S0 = X0*vk_2d(:,1);              % Compute S(i) using vk_2d
        F0 = X0*vk_2d(:,2);              % Compute F(i) using vk_2d
        C0 = (X0*vk_2d(:,3)).^(-1/gam);  % Compute C(i) using vk_2d 
        pie0 = ((1-(1-theta)*(S0/F0)^(1-epsil))/theta)^(1/(epsil-1));
                                         % Compute pie(i) from condition 
                                         % (35) in MM (2015)
        delta1 = ((1-theta)*((1-theta*pie0^(epsil-1))/(1-theta))^(epsil/(epsil-1))+theta*pie0^epsil/delta0)^(-1);
                                         % Compute delta(i) from condition 
                                         % (36) in MM (2015)
        Y0 = C0/(1-Gbar/exp(nuG0));      % Compute Y(i) from condition (38)
                                         % in MM (2015)
        L0 = Y0/exp(nua0)/delta1;        % Compute L(i) from condition (37) 
                                         % in MM (2015)
        Yn0 = (exp(nua0)^(1+vartheta)*(1-Gbar/exp(nuG0))^(-gam)/exp(nuL0))^(1/(vartheta+gam));
                                         %  Compute Yn(i) from condition (31) 
                                         % in MM (2015)
        R1 = piestar/betta*(R0*betta/piestar)^mu*((pie0/piestar)^phi_pie * (Y0/Yn0)^phi_y)^(1-mu)*exp(nuR0);   
                                         % Compute R(i) from conditions (27),
                                         % (39) in MM (2015)
        if zlb==1;R1=R1.*(R1>=1)+(R1<1); end
                                         % If ZLB is imposed, set R(i)=1 
                                         % if ZLB binds
        
        % 10.1.3 Next-period choices in grid point i
        %-------------------------------------------
        delta1_dupl = ones(n_nodes,1)*delta1; 
        R1_dupl = ones(n_nodes,1)*R1;
        % Duplicate "delta1" and "R1" n_nodes times to create a matrix with
        % n_nodes identical rows; n_G-by-n_nodes
                

        X1 = Ord_Polynomial_N([log(R1_dupl) log(delta1_dupl) nuR1 nua1 nuL1 nuu1 nuB1 nuG1],Degree);
        % Form complete polynomial of degree "Degree" on next-period state 
        % variables; n_nodes-by-npol_2d 
       
        S1 = X1*vk_2d(:,1);             % Compute next-period S using the 
                                        % corresponding vk_2d
        F1 = X1*vk_2d(:,2);             % Compute next-period F using the 
                                        % corresponding vk_2d
        C1 = (X1*vk_2d(:,3)).^(-1/gam); % Compute next-period C using the 
                                        % corresponding vk_2d 
        pie1 = ((1-(1-theta)*(S1./F1).^(1-epsil))/theta).^(1/(epsil-1));
                                        % Compute next-period pie using 
                                        % condition (35) in MM (2015)

        % 10.1.4. Evaluate conditional expectations in the Euler equations
        %------------------------------------------------------------------
        e(i,1) = (exp(nuu0)*exp(nuL0)*L0^vartheta*Y0/exp(nua0) + betta*theta*pie1.^epsil.*S1)'*weight_nodes; 
        e(i,2) = (exp(nuu0)*C0^(-gam)*Y0 + betta*theta*pie1.^(epsil-1).*F1)'*weight_nodes;
        e(i,3) = (betta*exp(nuB0)/exp(nuu0)*R1*exp(nuu1).*C1.^(-gam)./pie1)'*weight_nodes;

        % 10.1.5 Variables of the current iteration 
        %------------------------------------------
        S0_new_G(i,1) = S0(1,1);     
        F0_new_G(i,1) = F0(1,1);
        C0_new_G(i,1) = C0(1,1);

    end  
    
 % 10.2 Compute and update the coefficients of the policy functions 
 % ----------------------------------------------------------------
 vk_hat_2d = X0_G\e; % Compute the new coefficients of the 3 policy functions                          
 vk_2d(:,1) = damp(1)*vk_hat_2d(:,1) + (1-damp(1))*vk_2d(:,1); 
                     % Update the coefficients using damping
 vk_2d(:,2) = damp(2)*vk_hat_2d(:,2) + (1-damp(2))*vk_2d(:,2); 
                     % Update the coefficients using damping
 vk_2d(:,3) = damp(3)*vk_hat_2d(:,3) + (1-damp(3))*vk_2d(:,3); 
                     % Update the coefficients using damping
  
 % 10.3 Evaluate the percentage (unit-free) difference between the values on  
 % the grid from the previous and current iterations
 % -------------------------------------------------------------------------
 dif_EDS = mean(mean(abs(1-S0_new_G./S0_old_G)))/damp(1)+mean(mean(abs(1-F0_new_G./F0_old_G)))/damp(2)+mean(mean(abs(1-C0_new_G./C0_old_G)))/damp(3)
                   % The convergence criterion is adjusted to the damping 
                   % parameters   
                                          
 % 10.4. Store the values on the grid to be used on the subsequent iteration 
 %-------------------------------------------------------------------------
 S0_old_G = S0_new_G; 
 F0_old_G = F0_new_G; 
 C0_old_G = C0_new_G; 
end

% 11. Output of the iterative cycle
%----------------------------------
VK(:,:,2) = vk_2d;           % Fill in the matrix of coefficients obtained 
                             % in the main iterative cycle into VK

time_EDS = cputime - time0;  % Compute time used to solve the model

% 12. Simulate the EDS solution and evaluate its accuracy
%--------------------------------------------------------
 
% Generate the series for shocks using PER1
% -----------------------------------------
nua_test = PER1_test(4,2:end)'; 
nuL_test = PER1_test(7,2:end)';
nuR_test = PER1_test(8,2:end)';
nuG_test = PER1_test(9,2:end)';
nuB_test = PER1_test(10,2:end)';
nuu_test = PER1_test(11,2:end)';

% Initial values for the endogenous state variables
%--------------------------------------------------
R_init  = PER1_test(5,1);    % Nominal interest rate in the initial period
delta_init  = PER1_test(6,1);% Price dispersion in the initial period


[S_test F_test delta_test C_test Y_test Yn_test L_test R_test pie_test] = NK_EDS_simulation(vk_2d,nuR_test,nua_test,nuL_test,nuu_test,nuB_test,nuG_test,R_init,delta_init,gam,vartheta,epsil,betta,phi_y,phi_pie,mu,theta,piestar,Gbar,zlb,Degree);%end
                % Simulate the EDS solution

discard=200;    % The number of observations to discard
[Residuals_mean(5) Residuals_max(5) Residuals_max_E(1:9,5)] = NK_EDS_accuracy(nua_test,nuL_test,nuR_test,nuG_test,nuB_test,nuu_test,R_test,delta_test,L_test,Y_test,Yn_test,pie_test,S_test,F_test,C_test,rho_nua,rho_nuL,rho_nuR,rho_nuu,rho_nuB,rho_nuG,gam,vartheta,epsil,betta,phi_y,phi_pie,mu,theta,piestar,vcv,discard,vk_2d,Gbar,zlb,Degree);
                % Compute the mean and maximum residuals, "Residuals_mean" 
                % and "Residuals_max", across all points and all equilibrium
                %  conditions, as well as the maximum absolute  residuals, 
                % "Residuals_max_E", across all points, disaggregated  
                % by optimality conditions


Residuals_mean  % Display mean residuals; 1-by-5
Residuals_max   % Display maximum residuals; 1-by-5
Residuals_max_E % Display maximum residuals by equilibrium conditions; 9-by-5
    
save pie1_D3;   % Save the results
CPU = toc       % Save the elapsed time in CPU