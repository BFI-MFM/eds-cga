% The folder "OneCountryModel" contains MATLAB software that solves a 
% one-country model using an epsilon distinguishable set algorithm (EDS) 
% and cluster grid algorithm (CGA); it accompanies the article "Merging 
% simulation and projection approaches to solve high-dimensional problems 
% with an application to a new Keynesian model" by Lilia Maliar and Serguei 
% Maliar (2015), Quantitative Economics 6/1, pages 1–47 (henceforth, MM, 2015). 

% This software is based on that of Lilia  Maliar and Serguei Maliar for 
% solving models using the generalized stochastic simulation algorithm (GSSA) 
% method, as described in the paper "Numerically Stable and Accurate Stochastic
% Simulation Approaches for Solving Dynamic Economic Models" by Kenneth L. 
% Judd, Lilia Maliar and Serguei Maliar, (2011), Quantitative Economics 2/2, 
% 173–210 (henceforth, JMM, 2011). The modifications made are concerned with 
% the construction of the grid on which solutions are computed.
% 
% This version: March 19, 2015. First version: April 17, 2010.
%  
% 1. "Main_EDS_CGA_1.m"         is a main file for computing a solution to 
%                               the one-country model using the EDS and 
%                               CGA methods 
% 2. "Accuracy_Test_1.m"        computes residuals of the optimality 
%                               conditions of the one-country model on a   
%                               given set of points in the state space; 
%                               borrowed from JMM (2011)   
% 3. "Density.m"                estimates the density function from a 
%                               given set of points 
% 4. "Clusters.m"               constructs clusters from simulated series 
%                               and computes clusters' centers (to be used 
%                               as a grid) 
% 5. "EDS"                      constructs an epsilon distinguishable set 
%                               for a given set of data (to be used as a 
%                               grid)
% 6. "Ord_Polynomial_N.m"       constructs the sets of basis functions for                                 
%                               ordinary polynomials of the degrees from 
%                               one to five; borrowed  from JMM (2011)                          
% 7. "GH_Quadrature.m"          constructs integration nodes and weights for  
%                               the Gauss-Hermite rules with the number of  
%                               nodes in each dimension ranging from one to  
%                               ten; borrowed from JMM (2011)                     
% 8. "epsi10200.mat"            contains the series of the productivity  
%                               levels of length 10,200 that are used for 
%                               computing solutions and for evaluating  
%                               accuracy; borrowed from JMM (2011)
%
%
% -------------------------------------------------------------------------
% Copyright © 2015 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% The code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

clc;
clear all;

% 1. Choose the simulation length 
% ------------------------------- 
T     = 10000;   % Choose the simulation length for the solution procedure,
                 % T<=10,000                    
                 
% To solve models with T>10,000, one needs to simulate new series of the 
% productivity levels by enabling the code in paragraph 5  

% 2. Model's parameters
% ---------------------
gam     = 1;         % Utility-function parameter
alpha   = 0.36;      % Capital share in output
beta    = 0.99;      % Discount factor
delta   = 0.02;      % Depreciation rate 
rho     = 0.95;      % Persistence of the log of the productivity level
sigma   = 0.01;      % Standard deviation of shocks to the log of the 
                     % productivity level

% 3. Technology level
% -------------------
A       = (1/beta-(1-delta))/alpha;  
                     % Normalize steady state capital to one             

% 4. Initial condition
% --------------------
k(1:T+1,1) = 1;  % Initial condition for capital (is equal to steady state)
a(1:T,1)   = 1;  % Initial condition for the productivity level (is equal to 
                 % steady state)

% 5. Construct the productivity levels, a, for the solution procedure 
% -------------------------------------------------------------------
%epsi10200 = randn(10200,2);    % Generate a random draw of the productivity 
                                % shocks of length 10,000 
%save epsi10200 epsi10200;      % Save the series of the productivity shocks  
                                % into a file "epsi10000.mat" 
load epsi10200;                 % Load the previously saved series of the 
                                % productivity shocks of length 10,000 
epsi = epsi10200(1:T,1)*sigma;  % Compute the error terms in the process for 
                                % productivity level 
for t = 2:T; 
    a(t,1) = a(t-1,1)^rho*exp(epsi(t,1)); 
                                % Compute the next-period productivity levels
                                % using condition (7) in MM (2015)
end;

% _________________________________________________________________________
%                               
% Compute a first-degree polynomial solution using a GSSA method with one- 
% node Monte Carlo integration (this solution will be used as an initial 
% guess for the other cases) 
% _________________________________________________________________________
%
tic;                  % Start counting time needed to compute the solution
                           
% 6. The GSSA parameters  
% ---------------------
kdamp     = 0.01;     % Damping parameter for (fixed-point) iteration on 
                      % the coefficients of the capital policy function
dif_GSSA_1d = 1e+10;  % Set the initial difference between the series from
                      % two iterations in the convergence criterion (condition
                      % (10) in JMM, 2011) to a very large number

% 7. Initialize the first-degree capital policy function 
%-------------------------------------------------------                         
bk_1d    	= [0; 0.95; 0.05];  	          
                      % Vector of polynomial coefficients; 3-by-1
                      
% 8. Initialize the capital series
% --------------------------------
k_old = ones(T+1,1);   % Initialize the series of next-period capital; this
                       % series are used to check the convergence on the 
                       % subsequent iteration (initially, capital can take 
                       % any value); (T+1)-by-1

% 9. The main iterative cycle of GSSA
% -----------------------------------              
while dif_GSSA_1d > 1e-4*kdamp    % 1e-4*kdamp is a convergence parameter,
                                  % adjusted to the damping parameter; see 
                                  % JMM (2011) for a discussion
    
    % 9.1 Generate time series of capital
    % -----------------------------------
    for t = 1:T
      x(t,:) = [1 k(t,1) a(t,1)]; % The basis functions of the first-degree
                                  % polynomial at time t
      k(t+1,1) = x(t,:)*bk_1d;    % Compute next-period capital using 
                                  % polynomial coefficients bk_1d
    end;
   
    % 9.2 Compute time series of consumption, c, and Monte Carlo realizations 
    % of the right side of the Euler equation, y, defined in condition (8) 
    % in MM (2015)
    % -----------------------------------------------------------------------
    c    = A*k(1:T,1).^alpha.*a(1:T,1)  + (1-delta)*k(1:T,1)-k(2:T+1,1); 
                                 % T-by-1
    y    =  beta*c(2:T,1).^(-gam)./c(1:T-1,1).^(-gam).*(1-delta+alpha*A*k(2:T,1).^(alpha-1).*a(2:T,1)).*k(2:T,1);
                                 % (T-1)-by-1
                                 
    % 9.3 Evaluate the percentage (unit-free) difference between the series  
    % from the previous and current iterations
    % ---------------------------------------------------------------------
    dif_GSSA_1d = mean(abs(1-k./k_old))
                   % Compute a unit-free difference between the series from
                   % two iterations; see condition (10) in JMM (2011)
                                 
                                 
    % 9.4 Compute and update the coefficients of the capital policy function
    % ----------------------------------------------------------------------
    bk_hat_1d = inv(x(1:T-1,:)'*x(1:T-1,:))*x(1:T-1,:)'*y(1:T-1,:);
                                 % Compute new coefficients of the capital 
                                 % policy function using the OLS
    bk_1d = kdamp*bk_hat_1d + (1-kdamp)*bk_1d;  
                                 % Update the coefficients of the capital  
                                 % policy function using damping
                                     
    % 9.5 Store the capital series 
    %-----------------------------
    k_old = k;         % The stored capital series will be used for checking 
                       % the convergence on the subsequent iteration
    
end;

% 10. Time needed to compute the initial guess 
% --------------------------------------------
time_GSSA_1d     = toc; 

% _________________________________________________________________________                              
%
% Compute polynomial solutions of the degrees from one to D_max using either 
% EDS or CGA solution method 
% _________________________________________________________________________
%
tic;                  % Start counting time needed to compute the solution

% 11. The EDS_CGA parameters  
% -----------------------
kdamp = 0.05;            % Damping parameter for (fixed-point) iteration on 
                         % the coefficients of the capital policy function
dif_EDS_CGA_D = 1e+10;   % Set the initial difference between the series from
                         % two iterations in the convergence criterion
                         % to a very large number
                 
% 12. The matrices of the polynomial coefficients
% -----------------------------------------------                             
D_max  = 5;            % Maximum degree of a polynomial: the program computes
                       % polynomial solutions of the degrees from one to D_max;
                       % (D_max can be from 1 to 5) 
npol = [3 6 10 15 21]; % Number of coefficients in polynomials of the degrees                    
                       % from one to five
BK = zeros(npol(D_max),D_max); 
                       % Matrix of polynomial coefficients of the 
                       % capital policy function for the polynomial
                       % solutions of the degrees from one to D_max;
                       % npol(D_max)-by-D_max
BK(1:npol(1),1) = bk_1d;       
                       % Initial guess for the coefficients
                 
% 13. Choose an integration method for computing solutions  
% -------------------------------------------------------- 
IM  = 10;        % 0=a one-node Monte Carlo method (default);
                 % 1,2,..,10=Gauss-Hermite quadrature rules with 1,2,...,10 
                 % nodes, respectively
[n_nodes,epsi_nodes,weight_nodes] = GH_Quadrature(IM,1,sigma^2);
                 % n_nodes is the number of integration nodes, epsi_nodes  
                 % are integration nodes, and weight_nodes are integration 
                 % weights for Gauss-Hermite quadrature integration rule
                 % with IM nodes; for a one-dimensional integral, n_nodes  
                 % coincides with IM
a1 = a(1:T,:).^rho*exp(epsi_nodes');
                 % Compute the next-period productivity levels for each 
                 % integration node using condition (4) in JMM (2011); 
                 % T-by-n_nodes
                 
% 14. Construct a grid of representative points 
% ---------------------------------------------   
% A grid of representative points will be used to compute a solution to 
% the model. The grid is constructed from a given set of simulated points. 
% Two options are available: "Clusters.m" clusters simulated data into M 
% clusters and computes the centers of the clusters and "EDS.m" constructs 
% epsilon distinguishable subset of simulated points such that any two points 
% in this set are located at the distance at least epsilon. 

% 14.1 Estimate the density function and remove low-density points 
% ----------------------------------------------------------------
Data = [k(1:T,:) a];  % Construct data from the simulated series produced by GSSA   
[Density,PCn,di_min] = Density(Data,Data);
                      % Estimate the density function in all simulated points
Data_sort = sortrows([Density di_min Data],1);
                      % Sort the simulated points by the density function 
cutoff = round(T*0.01); 
                      % Cutoff level is 1% of points with the lowest density
Data1 = Data_sort(1+cutoff:end,3:end);
                      % Remove the low-density points with a given cutoff level

% 14.2 Construct a grid of representative points 
% ----------------------------------------------
M = 25;               % Target number of grid points (must be larger than the
                      % number of polynomial terms (coefficients)
GridMethod = 2;       % Choose a method for constructing the grid: 
                      % "1" is an EDS grid, and otherwise, it is a cluster 
                      % grid

% 14.2.1 Construct EDS with a target number of M grid points using bisection
% --------------------------------------------------------------------------  

if GridMethod == 1;  % GridMethod=1 means that we construct an EDS grid
    
    % Construct EDS sets with parameters "epsilon1" and "epsilon2" that 
    % contain no data points and all data points, respectively; these are 
    % the initial values that are necessary for bisection; see Algorithm 
    % M_bar of MM (2015)
    % -------------------------------------------------------------------
    PCn_sort = sortrows([Density di_min PCn],1); 
                   % Sort the principal components PCn by the density function
    r1 = min(sqrt(sum(PCn_sort(1+cutoff:end,3:end).^2,2)));
                   % The distance from the center to the closest point; a 
                   % ball with the radius r1 contains no data points
    r2 = max(sqrt(sum(PCn_sort(1+cutoff:end,3:end).^2,2)));
                   % The distance from the center to the furthest point; a 
                   % ball with the radius r2 contains all data points
    epsilon1 = r1/2/M;          
                   % The lower bound epsilon for bisection
    epsilon2 = 2*r2;            
                   % The upper bound epsilon for bisection
    
    % Construct EDS sets with "epsilon" that contains the target number of
    % M points (i.e., approximately)
    % ---------------------------------------------------------------------
    size_Grid_old = 0; 
                   % The number of points in the old grid, initially 0 
    Grid = zeros(M,2);  
                   % Grid; initially, it consists of M points   

    while abs(size_Grid_old-size(Grid,1))>0;  % Compute the difference in the 
                                              % number of points in the
                                              % new and old grids
        size_Grid_old = size(Grid,1);         % Store the number of points 
                                              % in the old grid
        epsilon = (epsilon1+epsilon2)/2;      % Compute epsilon for constructing
                                              % the new grid
        [Grid] = EDS(Data1,epsilon);          % construct the new grid
    if size(Grid,1)>M;                        % If the new grid has more points
                                              % than needed, ...
        epsilon1 = epsilon;                   % Increase the lower bound
    else                                      % Otherwise,  
        epsilon2 = epsilon;                   % Decrease the upper bound
    end
    end                                       % Proceed with bisection until
                                              % the grid converges, i.e., the 
                                              % number of grid points does
                                              % not change from one iteration 
                                              % to another
    
% 14.2.2 Construct a cluster grid with M grid points
% --------------------------------------------------
else
    [Grid] = Clusters(Data1,M); % Construct M clusters on "Data1" and compute
                                % the clusters' centers; 
end

% 14.3 The number of points in the grid
% -------------------------------------  

N_G = size(Grid,1);    % The actual number of points in the grid; for the
                       % EDS method, N_G can be different from M, and for
                       % the cluster grid method, N_G is exactly equal to M

                       
% 15. Compute the polynomial solutions of the degrees from one to D_max
% ---------------------------------------------------------------------
for D = 1:D_max
    tic
    
  % 15.1 Initialization of EDS_CGA
    % ---------------------------
    bk_D = BK(1:npol(D),1);% Compute the initial guess on the coefficients 
                           % equal to linear approximation
    
    X_G = Ord_Polynomial_N(Grid,D);   
                           % Construct the polynomial bases on the series 
                           % of state variables from the previously computed 
                           % time-series solution
    k_old = ones(N_G,1);   % Initialize the series of next-period capital; 
                           % this series is used to check the convergence 
                           % on the subsequent iteration (initially, 
                           % capital can take any value); N_G-by-1
    dif_EDS_CGA_D  = 1e+10;% Convergence criterion (initially is not satisfied)
    
    % 15.2 The main iterative cycle of EDS_CGA
    % -------------------------------------              
    while dif_EDS_CGA_D > 1e-4/10^D*kdamp;   
                           % 10^(-4-D)*kdamp is a convergence parameter, 
                           % adjusted to the polynomial degree D and the 
                           % damping parameter kdamp; see the discussion in 
                           % JMM (2011)
    

    % 15.2.2 Compute consumption series of all countries 
    %---------------------------------------------------
    k0  =  Grid(:,1);                        % Capital states  
    a0  =  Grid(:,2);                        % Productivity states  
    k1  =  X_G*bk_D;                         % Next-period capital 
    c   =  A*k0.^alpha.*a0 - k1+k0*(1-delta);% Consumption;

    % 15.2.3 Approximate the conditional expectations for t=1,...T-1 using 
    % the integration method chosen
    %----------------------------------------------------------------------
    %      
    % 15.2.3.2 Deterministic integration methods approximate the values of 
    % conditional expectations, y, in the Euler equation as a weighted average 
    % of the values of the integrand in the given nodes with the given weights 
    % ------------------------------------------------------------------------
    a1  =  a0.^rho*exp(epsi_nodes');   
    % Compute the next-period productivity levels for each integration       
    % node; n_nodes-by-1
                 
    k1_dupl = k1*ones(1,n_nodes);  
    % Duplicate k1 n_nodes times to create a matrix with 
    % n_nodes identical columns; N_G-by-n_nodes
            
    c_dupl = c(1:N_G,:)*ones(1,n_nodes);
    % Duplicate c n_nodes times to create a matrix with 
    % n_nodes identical columns; N_G-by-n_nodes 
    
    for i=1:n_nodes
        k2(:,i)  =  Ord_Polynomial_N([k1_dupl(:,i) a1(:,i)],D)*bk_D; 
    end
    % Compute capital of period t+2 (chosen at t+1) using the
    % capital policy functions; n_nodes-by-1 

    c1 =A*k1_dupl.^alpha.*a1 - k2+k1_dupl*(1-delta);
    % Compute next-period consumption in n_nodes; n_nodes-by-1
                       

    y  =  beta*c1.^(-gam)./c_dupl.^(-gam).*(1-delta+alpha*A*k1_dupl.^(alpha-1).*a1).*k1_dupl*weight_nodes;       
    % Compute condition expectation of the right side of the Euler equation
    % (8) in MM (2015); N_G-by-1
      
    % 15.2.4 Evaluate the percentage (unit-free) difference between the 
    % capital series from the previous and current iterations
    % -----------------------------------------------------------------
    D 
    dif_EDS_CGA_D = mean(mean(abs(1-k1./k_old)))
                % Compute a unit-free difference between the capital series 
                % from two iterations; see condition (10) in JMM (2011)

    % 15.2.5 Compute and update the coefficients of the capital policy 
    % functions 
    % ----------------------------------------------------------------
    bk_hat_D = X_G(1:N_G,:)\y(1:N_G,:); 
                              % Compute new coefficients of the capital 
                              % policy function using the chosen 
                              % approximation method
    bk_D = kdamp*bk_hat_D + (1-kdamp)*bk_D; 
                              % Update the coefficients of the capital  
                              % policy functions using damping 
                                     
    % 15.2.6 Store the capital series 
    %--------------------------------
    k_old = k1;       % The stored capital series will be used for checking 
                      % the convergence on the subsequent iteration

    end;

     % 15.2.7 The EDS_CGA output for the polynomial solution of degree D
     % --------------------------------------------------------------
     BK(1:npol(D),D) = bk_D;     % Store the coefficients of polynomial  
                                 % of degree D that approximates capital 
                                 % policy function 
                                      
     time_EDS_CGA(D) = toc;      % Time needed to compute the polynomial  
                                 % solution of degree D 
end                         

% 16. Accuracy test of the EDS_CGA solutions: residuals on a stochastic simulation 
% -----------------------------------------------------------------------------

   % 16.1 Specify a set of points on which the accuracy is evaluated
   %----------------------------------------------------------------
   T_test = 10200;           % Choose the simulation length for the test 
                             % on a stochastic simulation, T_test<=10,200
                             
   epsi_test = epsi10200(1:T_test,2); 
                             % Restrict the series of the productivity 
                             % levels for the test on a stochastic 
                             % simulation to the given T_test<=10,200  
                             % and N<=10                             
                               
   epsi_test = epsi_test(1:T_test,1)*sigma;     
                             % Compute the error terms in the process for 
                             % productivity level 
                             
   a_test(1,1) = 1;          % Initial condition for the productivity level
  
   for t = 2:T_test; 
       a_test(t,1) = a_test(t-1,1)^rho*exp(epsi_test(t,1));
                             % Compute the next-period productivity levels
                             % using the process for productivity (7) in MM
                             % (2015)
   end;
  
   k_test(1,1) = 1;          % Initial condition for capital (equal to steady state)
 
   % 16.2 Choose an integration method for evaluating the accuracy of solutions
   %---------------------------------------------------------------------------
   IM_test = 10;             % See paragraph 13 for the integration options

   % To implement the test on a stochastic simulation with T_test>10,200, one
   % needs to simulate new series of the productivity levels with larger T_test 
   % by enabling the code in paragraph 5

   % 16.3 Compute residuals on a stochastic simulation for the EDS_CGA polynomial 
   % solution of the degrees from one to D_max
   % -----------------------------------------------------------------------
   for D = 1:D_max
      
      % 16.3.1 Simulate the time series solution under the given capital-
      % policy-function coefficients, BK(:,:,D) with D=1,...,D_max 
      %------------------------------------------------------------------
      for t = 1:T_test
        X_test = Ord_Polynomial_N([k_test(t,1) a_test(t,1)],D);
        % The basis functions of a polynomial of degree D at time t
        k_test(t+1,1) = X_test*BK(1:npol(D),D);
        % Compute next-period capital using BK(1:npol(D),D)
      end
      
      % 16.3.2 Residuals across 10,200 points on a stochastic simulation
      %-----------------------------------------------------------------
      discard = 200;  % Discard the first 200 observations to remove the 
                      % effect of the initial conditions 
      [Residuals_mean(D) Residuals_max(D) time_test(D)] = Accuracy_Test_1(sigma,rho,beta,gam,alpha,delta,A,k_test(1:T_test,1),a_test(1:T_test,1),BK(1:npol(D),D),D,IM_test,discard);
      % Residuals_mean is the unit-free average absolute Euler equation error 
      % (in log10); 
      % Residuals_max is the unit-free maximum absolute Euler equation error 
      % (in log10)
  
   end;

% 17. Display the results for the polynomial solutions of the degrees from 
% one to five  
% ------------------------------------------------------------------------
disp(' '); disp('           EDS_CGA OUTPUT:'); disp(' '); 
disp('RUNNING TIME (in seconds):'); disp('');
disp('a) for computing the solution'); 
disp('      1        2         3         4        5');disp(time_EDS_CGA);
disp('b) for implementing the accuracy test'); 
disp('      1        2         3         4        5'); disp(time_test);
disp('APPROXIMATION ERRORS (log10):'); disp(''); 
disp('a) mean error in the Euler equation'); 
disp('      1        2         3         4        5');disp(Residuals_mean)
disp('b) max error in the Euler equation'); 
disp('      1        2         3         4        5');disp(Residuals_max)

%Table=[Residuals_mean' Residuals_max' time_EDS_CGA']
