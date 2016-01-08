% NK_PER12_simul.m is a routine for simulating a time-series perturbation   
% solution to the new Keynesian model considered in the article "Merging 
% simulation and projection approaches to solve high-dimensional problems 
% with an application to a new Keynesian model" by Lilia Maliar and Serguei 
% Maliar (2015), Quantitative Economics 6/1, pages 1–47 (henceforth, MM, 2015). 
% This routine does not use other routines: it computes the decision rules 
% using the matrices of coefficients.
% -------------------------------------------------------------------------
% Inputs:    "epsi" is the given vector of current-period shocks;
%            "SS", "del2", "A", "B", "C", "D", "E" are the matrices of 
%            coefficients from which the perturbation solution is constructed;
%            "sigma_nua", "sigma_nuL", "sigma_nuR", "sigma_nuu", "sigma_nuB",
%            and "sigma_nuG" are the parameters of the model;
%            "zlb" is a dummy parameter which is equal to 0 when ZLB is not
%            imposed and is equal to 1 when it is imposed
%
% Output:    "PER1", "PER2" are the first-order and second-order
%            perturbation solutions simulated for T periods
% -------------------------------------------------------------------------
% Copyright © 2015 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% The code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

function [PER1, PER2]  = NK_PER12_simulation(SS,del2,A,B,C,D,E,epsi,sigma_nua,sigma_nuL,sigma_nuR,sigma_nuu,sigma_nuB,sigma_nuG,T,zlb)

% The variables in the Dynare program "NK_Dynare.mod" are declared in the 
% following order: L Y Yn nua R delta nuL nuR nuG nuB nuu pie S F C; This 
% is also the decision rule (DR) order.

PER1 = SS*ones(1,T+1);  % Intialize PER1 by assuming steady state values; 15-by-(T+1)
PER2 = SS*ones(1,T+1);  % Intialize PER2 by assuming steady state values; 15-by-(T+1)

% Simulate the first-order perturbation solution
%-----------------------------------------------
for t=1:T
    DEV_all = PER1(:,t) - SS; % Compute a deviation of the vector of all 
                              % variables from the steady state, i.e., y(t)-SS
    DEV = DEV_all(4:11);      % Consider y(t)-SS for the state variables
    EPSI = [epsi(t,1)*sigma_nuR epsi(t,2)*sigma_nua epsi(t,3)*sigma_nuL epsi(t,4)*sigma_nuu epsi(t,5)*sigma_nuB epsi(t,6)*sigma_nuG]';
                              % Form a matrix of exogenous variables (shocks) 
                              % multiplied by the corresponding standard deviation
    PER1(:,t+1) = SS + A*DEV + B*EPSI; 
    % Compute PER1 (see a comment in Section 2 of NK_PER12_main)
 
    if zlb == 1;  PER1(5,t+1) = PER1(5,t+1)*(PER1(5,t+1)>=1)+1*(PER1(5,t+1)<1);end 
    % If the ZLB is imposed, the 5th variable, which is R, is set to 1
    % whenever R implied by the decision rule is smaller than 1 
end

% Simulate the second-order perturbation solution
%------------------------------------------------
for t = 1:T
    DEV_all = PER2(:,t)-SS;   % Compute a deviation of the vector of all 
                              % variables from the steady state, i.e., y(t)-SS
    DEV = DEV_all(4:11);      % Consider y(t)-SS for the state variables
    EPSI = [epsi(t,1)*sigma_nuR epsi(t,2)*sigma_nua epsi(t,3)*sigma_nuL epsi(t,4)*sigma_nuu epsi(t,5)*sigma_nuB epsi(t,6)*sigma_nuG]';
                              % Form a matrix of exogenous variables (shocks) 
                              % multiplied by the corresponding standard deviation
    DEV2 = kron(DEV,DEV);     % Compute a Kronecker product of of the vector 
                              % of state variables 
    EPSI2 = kron(EPSI,EPSI);  % Compute a Kronecker product of of the exogenous 
                              % variables
    DEVEPSI = kron(DEV,EPSI); % Compute a Kronecker product of the vector of state
                              % variables by the vector of exogenous variables 
    PER2(:,t+1) = SS+A*DEV+B*EPSI+0.5*del2+0.5*C*DEV2+0.5*D*EPSI2+E*DEVEPSI;
    % Compute PER2 (see a comment in Section 3 of NK_PER12_main)
    if zlb == 1;  PER2(5,t+1)=PER2(5,t+1)*(PER2(5,t+1)>=1)+1*(PER2(5,t+1)<1);end
    % If the ZLB is imposed, the 5th variable, which is R, is set to 1
    % whenever R implied by the decision rule is smaller than 1 
end
