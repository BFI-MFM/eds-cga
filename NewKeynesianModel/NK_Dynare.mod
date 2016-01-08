// Dynare routine that solves a new Keynesian model described in the article 
// "Merging simulation and projection approaches to solve high-dimensional 
// problems with an application to a new Keynesian model" by Lilia Maliar 
// and Serguei Maliar (2015), Quantitative Economics 6/1, pages 1–47 
// (henceforth, MM, 2015). 

// This version: March 6, 2015. First version: June 27, 2011.


// Block 1. Preamble

var L Y Yn  nua R delta nuL nuR nuG nuB nuu pie S F  C; //Endogenous variables

varexo eps_nuR eps_nua eps_nuL eps_nuu eps_nuB eps_nuG; // Exogenous variables

parameters gam vartheta epsil betta phi_y phi_pie mu rho_nua rho_nuL rho_nuR rho_nuu rho_nuB rho_nuG theta Gbar piestar sigma_nua sigma_nuL sigma_nuR sigma_nuu sigma_nuB sigma_nuG; // DECLARATION OF THE DEEP PARAMETERS.
                                                        // Parameters

// Load the parameters values
load parameterfile;

set_param_value('betta', betta);
set_param_value('piestar', piestar);
set_param_value('Gbar', Gbar);

set_param_value('gam', gam);
set_param_value('vartheta', vartheta);
set_param_value('epsil', epsil);
set_param_value('phi_y', phi_y);
set_param_value('phi_pie', phi_pie);
set_param_value('mu', mu);
set_param_value('theta', theta);

// autoregression coefficients in the AR(1) processes for shocks
set_param_value('rho_nua', rho_nua);
set_param_value('rho_nuu', rho_nuu);
set_param_value('rho_nuG', rho_nuG);
set_param_value('rho_nuL', rho_nuL);
set_param_value('rho_nuR', rho_nuR);
set_param_value('rho_nuB', rho_nuB);

// standard deviations of innovations in the AR(1) processes for shocks
set_param_value('sigma_nua', sigma_nua);
set_param_value('sigma_nuu', sigma_nuu);
set_param_value('sigma_nuG', sigma_nuG);
set_param_value('sigma_nuL', sigma_nuL);
set_param_value('sigma_nuR', sigma_nuR);
set_param_value('sigma_nuB', sigma_nuB);

// Block 2. Declaration of the model

model;

S = exp(nuu)*exp(nuL)*L^vartheta*Y/exp(nua) + betta*theta*pie(+1)^epsil*S(+1);
F = exp(nuu)*C^(-gam)*Y + betta*theta*pie(+1)^(epsil-1)*F(+1);
S = ((1-theta*pie^(epsil-1))/(1-theta))^(1/(1-epsil))*F;
delta = ((1-theta)*((1-theta*pie^(epsil-1))/(1-theta))^(epsil/(epsil-1)) + theta*pie^epsil/delta(-1))^(-1);
C = (betta*exp(nuB)*exp(nuu(+1))/exp(nuu)*C(+1)^(-gam)*R/pie(+1))^(-1/gam);
Y = exp(nua)*L*delta;
C = (1-Gbar/exp(nuG))*Y;
Yn = (exp(nua)^(1+vartheta)*(1-Gbar/exp(nuG))^(-gam)/exp(nuL))^(1/(vartheta+gam));
R = piestar/betta*(R(-1)*betta/piestar)^mu*((pie/piestar)^phi_pie * (Y/Yn)^phi_y)^(1-mu)*exp(nuR);     //Taylor rule

nua = rho_nua*nua(-1) + eps_nua;
nuL = rho_nuL*nuL(-1) + eps_nuL;
nuR = rho_nuR*nuR(-1) + eps_nuR;
nuu = rho_nuu*nuu(-1) + eps_nuu;
nuB = rho_nuB*nuB(-1) + eps_nuB;
nuG = rho_nuG*nuG(-1) + eps_nuG;
end;

// Block 3. Solving the model
// Numerical initial conditions for the computation of the deterministic steady state 

initval;
Yn = exp(Gbar)^(gam/(vartheta+gam));
Y = Yn;
pie = 1;  // also try pie=piestar
delta = 1;
L = Y/delta;
C = (1-Gbar)*Y;
F = C^(-gam)*Y/(1-betta*theta*pie^(epsil-1));
//S = ((1-theta*pie^(epsil-1)) / (1-theta)  ) ^ (1/(1-epsil)) * F;
S =   L^vartheta*Y/(1-betta*theta*pie^epsil);
R = pie/betta;

nua = 0;
nuL = 0;
nuR = 0;
nuu = 0;
nuB = 0;
nuG = 0;
eps_nua = 0;
eps_nuL = 0;
eps_nuR = 0;
eps_nuu = 0;
eps_nuB = 0;
eps_nuG = 0;
end;

// Specify shocks
shocks;
var eps_nua; stderr sigma_nua;
var eps_nuL; stderr sigma_nuL;
var eps_nuR; stderr sigma_nuR;
var eps_nuu; stderr sigma_nuu;
var eps_nuB; stderr sigma_nuB;
var eps_nuG; stderr sigma_nuG;
end;

// The model is solved and simulated with the "stoch_simul" command
stoch_simul(order=2, nograph, noprint, nomoments, irf=0);

