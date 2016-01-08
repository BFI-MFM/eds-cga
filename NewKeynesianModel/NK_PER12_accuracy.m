% NK_PER12_accuracy.m is a routine for evaluating accuracy of the perturbation     
% solutions to the new Keynesian model considerd in the article "Merging 
% simulation and projection approaches to solve high-dimensional problems 
% with an application to a new Keynesian model" by Lilia Maliar and Serguei 
% Maliar (2015), Quantitative Economics 6/1, pages 1–47 (henceforth, MM, 2015). 
% It computes residuals in the optimality conditions on a given set of points 
% in the state space.
% -------------------------------------------------------------------------
% Inputs:    "PER" is a perturbation solution in the given set of points on  
%            which the accuracy is tested; 
%            "SS", "del2", "A", "B", "C", "D", "E" are the matrices of 
%            coefficients from which the perturbation solution is constructed;
%            "mu", "gam", "vartheta", "beta", "A", "tau", "rho" and "vcv"
%            are the parameters of the model;
%            "discard" is the number of data points to discard 
%
% Outputs:   "Residuals_mean" and "Residuals_max" are, respectively, the mean 
%            and maximum absolute residuals across all points and all equi-
%            librium conditions; 
%            "Residuals_max_E" is the maximum absolute residuals across all 
%            points, disaggregated  by optimality conditions; 
%            "Residuals" are absolute residuals disaggregated by the equi-
%            librium conditions
% -------------------------------------------------------------------------
% Copyright © 2015 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% The code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

function [Residuals_mean Residuals_max Residuals_max_E Residuals] = NK_PER12_accuracy(SS,del2,A,B,C,D,E,PER,order,gam,vartheta,epsil,beta,phi_y,phi_pie,mu,theta,Gbar,piestar,vcv,discard,zlb)

tic                   % Start counting time for running the test        

T = size(PER,2);      % Infer the number of points on which accuracy is 
                      % evaluated 
                      
% Integration method for evaluating accuracy 
% ------------------------------------------
[n_nodes,epsi_nodes,weight_nodes] = Monomials_2(6,vcv);
                             % Monomial integration rule with 2N^2+1 nodes;
                             % there are 6 shocks
                             
% Given the simulated time-series perturbation solution, call the variables 
% by names
%--------------------------------------------------------------------------
% The variables in the Dynare program "NK_Dynare.mod" are declared in the 
% following order: L Y Yn nua R delta nuL nuR nuG nuB nuu pie S F C; This is 
% also the decision rule (DR) order

PER_nodes = SS*zeros(1,n_nodes);% Intialize PER_nodes (perturbation solution)
                                % in all future nodes by assuming steady state 
                                % values; 15-by-n_nodes

for t = 2:T;              % For each given point, ...     
        t
        nua0 = PER(4,t);  % nua is the 4th variable in the declaration order
        nuL0 = PER(7,t);  % nuL0 is the 7th variable in the declaration order
        nuR0 = PER(8,t);  % nuR is the 8th variable in the declaration order
        nuG0 = PER(9,t);  % nuG is the 9th variable in the declaration order
        nuB0 = PER(10,t); % nuB is the 10th variable in the declaration order
        nuu0 = PER(11,t); % nuu is the 11th variable in the declaration order

        R0 = PER(5,t-1);  % R(t-1) is the 5th variable in the declaration order
                          % but we taken one period before (this is a state 
                          % of period t)
        delta0 = PER(6,t-1);% delta(t-1) is the 6th variable in the declaration 
                          % order taken one period before (this is a state 
                          % of period t)
        R1 = PER(5,t);    % R(t) is the 5th variable in the declaration order 
        delta1 = PER(6,t);% delta(t) is the 6th variable in the declaration 
                          % order 
        L0 = PER(1,t);    % L is the 1st variable in the declaration order
        Y0 = PER(2,t);    % Y is the 2d variable in the declaration order
        Yn0 = PER(3,t);   % Yn is the 3d variable in the declaration order
        pie0 = PER(12,t); % pie is the 12th variable in the declaration order
        S0 = PER(13,t);   % S is the 13th variable in the declaration order
        F0 = PER(14,t);   % F is the 14th variable in the declaration order
        C0 = PER(15,t);   % C is the 15th variable in the declaration order
       
% 3. Compute the values of all the variables in all future nodes using the 
% matrices of coefficients SS, de l2, A, B, C, D, E (perturbation solution)
% ------------------------------------------------------------------------      
if order == 1              % Consider the first-order perturbation solution
    DEV_all = PER(:,t)-SS; % Compute a deviation of the vector of all 
                           % variables from the steady state, i.e., y(t)-SS
    DEV = DEV_all(4:11);   % Consider y(t)-SS for the state variables
    for i = 1:n_nodes;
    EPSI = [epsi_nodes(i,1) epsi_nodes(i,2) epsi_nodes(i,3) epsi_nodes(i,4) epsi_nodes(i,5) epsi_nodes(i,6)]';
    % Form a matrix of shocks
    PER_nodes(:,i) = SS+A*DEV+B*EPSI;
    % Compute the first-order perturbation solution in the future nodes 
    % (see a comment in Section 3.1 of "Main_NK_Degree2")
    end
elseif order == 2          % Consider the second-order perturbation solution
    DEV_all = PER(:,t)-SS; % Compute a deviation of the vector of all 
                           % variables from the steady state, i.e., y(t)-SS
    DEV = DEV_all(4:11);   % Consider y(t)-SS for the state variables
   for i = 1:n_nodes 
    EPSI = [epsi_nodes(i,1) epsi_nodes(i,2) epsi_nodes(i,3) epsi_nodes(i,4) epsi_nodes(i,5) epsi_nodes(i,6)]';
    % Form a matrix of shocks
    DEV2 = kron(DEV,DEV);    % Compute a Kronecker product of the vector 
                             % of state variables 
    EPSI2 = kron(EPSI,EPSI); % Compute a Kronecker product of the exogenous 
                             % variables
    DEVEPSI = kron(DEV,EPSI);% Compute a Kronecker product of the vector of state
                             % variables by the vector of exogenous variables 
    PER_nodes(:,i) = SS+A*DEV+B*EPSI+0.5*del2+0.5*C*DEV2+0.5*D*EPSI2+E*DEVEPSI;
    % Compute the second-order perturbation solution in the future nodes  
    % (see a comment in Section 3.2 of "Main_NK_Degree2")
   end
end

   % Call the future variables by names
   %-----------------------------------
   pie1 = PER_nodes(12,:);
   S1   = PER_nodes(13,:);
   F1   = PER_nodes(14,:);
   C1   = PER_nodes(15,:);
   nuu1 = PER_nodes(11,:);

    % Compute residuals for each of the 9 equilibrium conditions
    %-----------------------------------------------------------
    Residuals(t-1,1) = 1-(exp(nuu0)*exp(nuL0)*L0^vartheta*Y0/exp(nua0) + beta*theta*pie1.^epsil.*S1)*weight_nodes/S0; 
    Residuals(t-1,2) = 1-(exp(nuu0)*C0^(-gam)*Y0 + beta*theta*pie1.^(epsil-1).*F1)*weight_nodes/F0;
    Residuals(t-1,3) = 1-(beta*exp(nuB0)/exp(nuu0)*R1*exp(nuu1).*C1.^(-gam)./pie1)*weight_nodes/C0^(-gam);
    Residuals(t-1,4) = 1-((1-theta*pie0^(epsil-1))/(1-theta))^(1/(1-epsil))*F0/S0;
    Residuals(t-1,5) = 1-((1-theta)*((1-theta*pie0^(epsil-1))/(1-theta))^(epsil/(epsil-1)) + theta*pie0^epsil/delta0)^(-1)/delta1;
    Residuals(t-1,6) = 1-exp(nua0)*L0*delta1/Y0;
    Residuals(t-1,7) = 1-(1-Gbar/exp(nuG0))*Y0/C0;
    Residuals(t-1,8) = 1-(exp(nua0)^(1+vartheta)*(1-Gbar/exp(nuG0))^(-gam)/exp(nuL0))^(1/(vartheta+gam))/Yn0;
    Residuals(t-1,9) = 1-piestar/beta*(R0*beta/piestar)^mu*((pie0/piestar)^phi_pie * (Y0/Yn0)^phi_y)^(1-mu)*exp(nuR0)/R1;   % Taylor rule
    if zlb==1; Residuals(t-1,9) = Residuals(t-1,9)*(R1>1);end
    % If the ZLB is imposed and R>1, the residuals in the Taylor rule (the 
    % 9th equation) are zero
end

% Residuals across all the equilibrium conditions and all test points
%--------------------------------------------------------------------
 Residuals_mean = log10(mean(mean(abs(Residuals(1+discard:end,:))))); 
 % Mean absolute residuals computed after discarding the first 
 % "discard" observations
 
 Residuals_max = log10(max(max(abs(Residuals(1+discard:end,:)))));    
 % Maximum absolute residuals computed after discarding the first 
 % "discard" observations
 
 % Residuals disaggregated by the eqiulibrium conditions
 %------------------------------------------------------
 Residuals_max_E = log10(max(abs(Residuals(1+discard:end,:))))';    
 % Maximum absolute residuals across all test points for each of the 9 
 % equilibrium conditions computed after discarding the first "discard" 
 % observations; 
 
time_test = toc;     % Time needed to run the test