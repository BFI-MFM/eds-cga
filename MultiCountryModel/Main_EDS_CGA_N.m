% The folder "Main_EDS-CGA_N.m contains MATLAB software that solves a multi-
% country neoclassical growth model using an epsilon distinguishable set 
% algorithm (EDS) and cluster grid algorithm (CGA), as described in the 
% article "Merging simulation and projection approaches to solve high-
% dimensional problems with an application to a new Keynesian model" by 
% Lilia Maliar and Serguei Maliar (2015), Quantitative Economics 6/1, 
% pages 1–47 (henceforth, MM, 2015). 
 
% This software is based on that of Lilia  Maliar and Serguei Maliar for 
% solving the multi-country model using the generalized stochastic simulation 
% algorithm (GSSA) method, as described in the paper "Numerically Stable and 
% Accurate Stochastic Simulation Approaches for Solving Dynamic Economic 
% Models" by Kenneth L. Judd, Lilia Maliar and Serguei Maliar, (2011), 
% Quantitative Economics 2/2, 173–210 (henceforth, JMM, 2011). The modifi-
% cations made are concerned with the construction of the grid on which 
% solutions are computed. 
% 
% This version: March 19, 2015. First version: April 17, 2010.
% 
% ------------------------------------------------------------------------
% The software uses the following files: 
% ------------------------------------------------------------------------
% 1. "Main_EDS_CGA_N.m"          is a main file for computing a solution to 
%                                the multi-country model using the EDS and 
%                                CGA methods
% 2. "Accuracy_Test_N.m"         computes residuals of the optimality 
%                                conditions of the multi-country model on a   
%                                given set of points in the state space; 
%                                borrowed from JMM (2011)  
% 3. "Density.m"                 estimates the density function from a 
%                                given set of points 
% 4. "Clusters.m"                constructs clusters from simulated series 
%                                and computes clusters' centers (to be used 
%                                as a grid) 
% 5. "EDS"                       constructs an epsilon distinguishable set 
%                                for a given set of data (to be used as a 
%                                grid)
% 6. "Ord_Polynomial_N.m"        constructs the sets of basis functions for 
%                                ordinary polynomials of the degrees from 
%                                one to five, for the N-country model;  
%                                borrowed from JMM (2011)
% 7. "Productivity.m"            generates random draws of the productivity  
%                                shocks and simulates the corresponding   
%                                series of the productivity levels; borrowed  
%                                from JMM (2011)
% 8. "Monomials_1.m"             constructs integration nodes and weights 
%                                for an N-dimensional monomial (non-product) 
%                                integration rule with 2N nodes; borrowed  
%                                from JMM (2011) 
% 9. "Monomials_2.m"             constructs integration nodes and weights for 
%                                an N-dimensional monomial (non-product) 
%                                integration rule with 2N^2+1 nodes; borrowed 
%                                from JMM (2011)
% 10. "GH_Quadrature.m"          constructs integration nodes and weights for  
%                                the Gauss-Hermite rules with the number of  
%                                nodes in each dimension ranging from one to  
%                                ten; borrowed from JMM (2011)                     
% 11. "aT20200N10.mat"           contains the series of the productivity  
%                                levels of length 20,200 for 10 countries that 
%                                are used for computing solutions and for 
%                                evaluating accuracy; borrowed from JMM (2011)
% -------------------------------------------------------------------------
% Copyright © 2015 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% The code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

clc;
clear all;

% 1. Choose the number of countries and simulation length
% -------------------------------------------------------
N      = 2;           % Choose the number of countries, 1<=N<=10 (note that the 
                      % code also works for the one-country case, N=1)
D_max  = 5;           % Maximum degree of a polynomial: the program computes
                      % polynomial solutions of the degrees from one to D_max;
                      % (D_max can be from 1 to 5 but for D_max>2, smaller 
                      % dumping parameter kdamp may be needed for convergence) 
npol_2d = size(Ord_Polynomial_N(ones(1,2*N),D_max),2);
                      % Number of terms in a complete given degree polynomial 
                      % for the case of N countries
T     = 10000;        % Choose the simulation length for the solution procedure

% To solve models with N>10 or T>10,000, one needs to simulate new series
% of the productivity levels by enabling the code in paragraph 5  

% 2. Model's parameters
% ---------------------
gam     = 1;      % Utility-function parameter
alpha   = 0.36;   % Capital share in output
beta    = 0.99;   % Discount factor
delta   = 0.025;  % Depreciation rate 
rho     = 0.95;   % Persistence of the log of the productivity level
sigma   = 0.01;   % Standard deviation of shocks to the log of the 
                  % productivity level
vcv = sigma^2*(eye(N)+ones(N,N)); 
                  % Variance-covariance matrix of the countries' productivity 
                  % shocks in which diagonal terms are equal to 2*sigma^2   
                  % and in which off-diagonal terms are equal to sigma^2; 
                  % this vcv follows from the assumption that a country's 
                  % shock has both common-for-all-countries and country-
                  % specific components; N-by-N

                                     
% 3. The normalizing constant, A, and welfare weight, tau
% -------------------------------------------------------
A       = (1-beta+beta*delta)/alpha/beta;  % The normalizing constant in output  
tau     = 1;                               % The welfare weight of country j 

% The above normalization ensures that steady state of capital of all 
% countries is equal to one 

% 4. Initial condition
% --------------------
k(1,1:N) = 1;  % Initial condition for capital (is equal to steady state)
a(1,1:N) = 1;  % Initial condition for the productivity level (is equal to 
               % steady state)

% 5. Construct the productivity levels, a 
% ---------------------------------------
% a20200 = Productivity(T,N,a(1,1:N),sigma,rho);
                               % Generate a random draw of the productivity 
                               % shocks and simulate the corresponding  
                               % series of the productivity levels of length   
                               % T periods for N countries 
% save aT20200N10 a20200;       % Save the series of the productivity levels  
                               % into a file "aT20200N10.mat" 
load aT20200N10;               % Load the previously saved series of the 
                               % productivity levels of length 20,200 for 
                               % 10 countries (the first 10,000 observations
                               % are used for finding a solution, and the 
                               % remaining 10,200 observations are used for
                               % evaluating accuracy)
a = a20200(1:T,1:N);           % Restrict the series of the productivity 
                               % levels for the solution procedure to the 
                               % given T<=10,000 and N<=10
% _________________________________________________________________________
%                               
% Compute a first-degree polynomial solution using the GSSA method with one-
% node Monte Carlo integration method (this solution will be used as an initial 
% guess for the EDS and CGA methods) 
% _________________________________________________________________________
%
tic;                  % Start counting time needed to compute the solution
                            
% 6. The GSSA parameters  
% ---------------------
kdamp     = 0.1;     % Damping parameter for (fixed-point) iteration on 
                     % the coefficients of the capital policy functions
dif_1d = 1e+10;      % Set the initial difference between the series from
                     % two iterations in the convergence criterion (condition
                     % (10) in JMM, 2011) to a very large number
% To achieve convergence under N>10, one may need to modify the values of 
% the damping parameter kdamp or refine the initial guess 

% 7. Initialize the first-degree capital policy functions of N countries 
%-----------------------------------------------------------------------                          
bk_1d  = [zeros(1,N); diag(ones(1,N)*0.9);diag(ones(1,N)*0.1)]; 
% Matrix of polynomial coefficients of size (1+2N)-by-N: for each country  
% (each column), 1+2N rows correspond to a constant, N coefficients on the
% capital stocks, k(t,1),...,k(t,N), and N coefficients on the productivity 
% levels, a(t,1),...,a(t,N)

% As an initial guess, assume that a country's j capital depends only on 
% its own capital and productivity level as k(t+1,j)=0.9*k(t,j)+0.1*a(t,j); 
% (note that in the steady state, we have k(t+1,j)=0.9*k(t,j)+0.1*a(t,j)=1)

% Note that diag(ones(1,N)*q) delivers an N-by-N matrix with  diagonal 
% entries equal to q. 

% 8. Initialize the capital series
% --------------------------------
k_old = ones(T+1,N);   % Initialize the series of next-period capital of N
                       % countries; these series are used to check the
                       % convergence on the subsequent iteration (initially, 
                       % capital can take any value); (T+1)-by-N

x = zeros(T,2*N+1);
                       
% 9. The main iterative cycle of GSSA
% -----------------------------------              
while dif_1d > 1e-4*kdamp;        % 10^4*kdamp is a convergence parameter,
                                  % adjusted to the damping parameter; see 
                                  % JMM (2011) for a discussion

% 9.1 Generate time series of capital
% -----------------------------------
for t = 1:T 
    x(t,:) = [1 k(t,:) a(t,:)];   % The basis functions of the first-degree 
                                  % polynomial of at time t
    k(t+1,:) = x(t,:)*bk_1d;      % Compute next-period capital using 
                                  % polynomial coefficients bk_1d
end

% 9.2 Compute consumption series 
%-------------------------------
C = (A*k(1:T,:).^alpha.*a(1:T,:) - k(2:T+1,:)+k(1:T,:)*(1-delta))*ones(N,1);
% Aggregate consumption is computed by summing up individual consumption, 
% which in turn, is found from the individual budget constraints; T-by-1

c = C*ones(1,N)/N;                % Individual consumption is the same for
                                  % all countries; T-by-N 

% 9.3 Evaluate the percentage (unit-free) difference between the series  
% from the previous and current iterations
% ---------------------------------------------------------------------
dif_1d = mean(mean(abs(1-k./k_old)))
                % Compute a unit-free difference between the capital series 
                % from two iterations; see condition (10) in JMM (2011)
                   
% 9.4 Monte Carlo realizations of the right side of the Euler equation, Y, 
% in condition (12) in MM (2015)
%-------------------------------------------------------------------------
for j = 1:N
   Y(1:T-1,j) = beta*c(2:T,j).^(-gam)./c(1:T-1,j).^(-gam).*(1-delta+alpha*A*k(2:T,j).^(alpha-1).*a(2:T,j)).*k(2:T,j);
end  % (T-1)-by-N

% 9.5 Compute and update the coefficients of the capital policy functions 
% -----------------------------------------------------------------------
bk_hat_1d  = x(1:T-1,:)\Y(1:T-1,:);    
                              % Compute new coefficients of the capital 
                              % policy functions using the backslash
                              % operator
bk_1d = kdamp*bk_hat_1d + (1-kdamp)*bk_1d; 
                              % Update the coefficients of the capital  
                              % policy functions using damping 
                                     
% 9.6 Store the capital series 
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

tic;                  % Start counting time needed to compute the solution
                            
% 11. The EDS and CGA parameters  
% ----------------------------
kdamp     = 0.1;      % Damping parameter for (fixed-point) iteration on 
                      % the coefficients of the capital policy functions
dif_EDS_D  = 1e+10;   % Set the initial difference between the series from
                      % two iterations in the convergence criterion to a 
                      % very large number

% To achieve convergence under N>10, one may need to modify the values of 
% the damping parameter kdamp or refine the initial guess 

% 12. The matrix of the polynomial coefficients
% ---------------------------------------------                             
for D = 1:D_max       % For the polynomial degrees from one to D_max
    npol(D) = size(Ord_Polynomial_N([k(1,:) a(1,:)],D),2);      
end
                      % Construct polynomial bases and compute their number
                      % (this is needed for finding the number of the polynomial
                      % coefficients in the policy functions)
                      
BK = zeros(npol(D_max),N,D_max);   % Matrix of polynomial coefficients of the 
                                   % capital policy functions for the polynomial
                                   % solutions of the degrees from one to D_max;
                                   % npol(D_max)-by-N-by-D_max
BK(1:npol(1),1:N,1) = bk_1d;       % Initial guess for the coefficients
 
% 13. Choose an integration method 
% --------------------------------                             
IM    = 12;      % 0=a one-node Monte Carlo method(default);
                 % 1,2,..,10=Gauss-Hermite quadrature rules with 1,2,...,10 
                 % nodes in each dimension, respectively;
                 % 11=Monomial rule with 2N nodes;
                 % 12=Monomial rule with 2N^2+1 nodes
if (IM>=1)&&(IM<=10)
    [n_nodes,epsi_nodes,weight_nodes] = GH_Quadrature(IM,N,vcv);
                             % Compute the number of integration nodes, 
                             % n_nodes, integration nodes, epsi_nodes, and 
                             % integration weights, weight_nodes, for Gauss-
                             % Hermite quadrature integration rule with IM 
                             % nodes in each dimension
elseif IM == 11
    [n_nodes,epsi_nodes,weight_nodes] = Monomials_1(N,vcv);
                             % Monomial integration rule with 2N nodes
elseif IM == 12
    [n_nodes,epsi_nodes,weight_nodes] = Monomials_2(N,vcv);
                             % Monomial integration rule with 2N^2+1 nodes
end

% Under the one-node Gauss-Hermite quadrature rule, the conditional
% expectation (integral) is approximated by the value of the integrand 
% evaluated in one integration node in which the next-period productivity 
% shock is zero, i.e., the next-period productivity level is 
% a(t+1,:)=a(t,:).^rho*exp(0)=a(t,:).^rho

    
% 14. Construct a grid of representative points 
% ---------------------------------------------   
% A grid of representative points will be used to compute a solution to 
% the model. The grid is constructed from a given set of simulated points. 
% Two options are available: "Clusters.m" clusters simulated data into M 
% clusters and computes the centers of the clusters, and "EDS.m" constructs 
% epsilon distinguishable subset of simulated points such that any two points 
% in this set are located at the distance at least epsilon. 

% 14.1 Estimate the density function and remove low-density points 
% ----------------------------------------------------------------
Data = [k(1:T,:) a];  % Construct data from the simulated series produced 
                      % by GSSA   
[Density,PCn,di_min] = Density(Data,Data);
                      % Estimate the density function in all simulated points
Data_sort = sortrows([Density di_min Data],1);
                      % Sort the simulated points by the density function 
cutoff = round(T*0.01); 
                      % Cutoff level is 1% of points with the lowest density
Data1=Data_sort(1+cutoff:end,3:end);
                      % Remove the low-density points with a given cutoff level

% 14.2 Construct a grid of representative points 
% ----------------------------------------------
M = 400;              % Target the number of grid points (must be larger 
                      % than the number of polynomial terms (coefficients)
GridMethod = 1;       % Choose a method for constructing the grid: 
                      % "1" is an EDS grid, and otherwise, it is a cluster 
                      % grid. 

% 14.2.1 Construct EDS with a target number of M grid points using bisection
% --------------------------------------------------------------------------  

if GridMethod==1;  % GridMethod=1 means that we construct an EDS grid
    
    % Construct EDS sets with parameters "epsilon1" and "epsilon2" that 
    % contain no data points and all data points, respectively; these are 
    % the initial values that are necessary for bisection; see Algorithm 
    % M_bar of MM (2015)
    % ---------------------------------------------------------------------
    PCn_sort = sortrows([Density di_min PCn],1); 
                   % Sort the principal components PCn by density function
    r1 = min(sqrt(sum(PCn_sort(1+cutoff:end,3:end).^2,2)));
                   % The distance from the center to the closest point; a  
                   % ball with the radius r1 contains no data points
    r2 = max(sqrt(sum(PCn_sort(1+cutoff:end,3:end).^2,2)));
                   % The distance from the center to the furthest point; a  
                   % ball with the radius r2 contains all data points
    epsilon1 = r1/2/M^(1/N);          
                   % The lower bound epsilon for bisection
    epsilon2 = 2*r2; 
                   % The upper bound epsilon for bisection
    
    % Construct an EDS sets with "epsilon" that contains the target number 
    % of M points (i.e., approximately)
    % ---------------------------------------------------------------------
    size_Grid_old = 0;    
                   % The number of points in the old grid, initially 0 
    Grid = zeros(M,N*2);  
                   % Grid; initially, it consists of M points   

    while abs(size_Grid_old-size(Grid,1))>0;  % Compute the difference in the 
                                              % number of points in the
                                              % new and old grids
        size_Grid_old = size(Grid,1);         % Store the number of points 
                                              % in the old grid
        epsilon = (epsilon1 + epsilon2)/2;    % Compute epsilon for constructing
                                              % the new grid
        [Grid] = EDS(Data1,epsilon);          % construct the new grid
    if size(Grid,1)>M;                        % If the new grid has more points
                                              % than needed, ...
        epsilon1 = epsilon;                   % Then, increase the lower bound
    else                                      % Otherwise, ...
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
    [Grid] = Clusters(Data1,M); % Construct M clusters on "Data1" and 
                                % compute the clusters' centers
end

% 14.3 The number of points in the grid
% -------------------------------------  

N_G = size(Grid,1);    % The actual number of points in the grid; for the
                       % EDS method, N_G can be different from M, and for
                       % the cluster grid method, N_G is exactly equal M

% 15. Compute the polynomial solutions of the degrees from one to D_max
% ---------------------------------------------------------------------
for D = 1:D_max
    tic
    
    % 15.1 Initialization of EDS and CGA
    % ----------------------------------
    bk_D = BK(1:npol(D),1:N,1);% Compute the initial guess on the coefficients 
                               % equal to linear approximation
                               
    X_G = Ord_Polynomial_N(Grid,D);   
                           % Construct the polynomial bases on the series 
                           % of state variables from the previously computed 
                           % time-series solution
    k_old = ones(N_G,N);   % Initialize the series of next-period capital of N
                           % countries; these series are used to check the
                           % convergence on the subsequent iteration (initially, 
                           % capital can take any value); N_G-by-N
    dif_EDS_CGA_D  = 1e+10;% Convergence criterion (initially is not satisfied)
    
    % 15.2 The main iterative cycle of EDS and CGA
    % --------------------------------------------              
    while dif_EDS_CGA_D > 1e-4/10^D*kdamp;   
                           % 10^(-4-D)*kdamp is a convergence parameter, 
                           % adjusted to the polynomial degree D and the 
                           % damping parameter kdamp; see the discussion in 
                           % JMM (2011)
    

    % 15.2.2 Compute consumption series of all countries 
    %---------------------------------------------------
    k0  =  Grid(:,1:N);                     % N current capital stocks  
    a0  =  Grid(:,N+1:2*N);                 % N current productivity levels  
    k1  =  X_G*bk_D;                        % N next-period capital stocks 
    C = (A*k0.^alpha.*a0 - k1+k0*(1-delta))*ones(N,1);
    % Aggregate consumption is computed by summing up individual consumption, 
    % which in turn, is found from the individual budget constraints; T-by-1


    c = C*ones(1,N)/N;            % Individual consumption is the same for
                                  % all countries; T-by-N 
 
    % 15.2.3 Approximate the conditional expectations for t=1,...T-1 using 
    % the integration method chosen
    %----------------------------------------------------------------------
    %      
    % 15.2.3.2 Deterministic integration methods approximate the values of 
    % conditional expectations, Y, in the Euler equation as a weighted average 
    % of the values of the integrand in the given nodes with the given weights 
    % ------------------------------------------------------------------------
        Y = zeros(N_G,N); % Allocate memory for the variable Y
        for i = 1:n_nodes         
            a1  =  Grid(:,N+1:2*N).^rho.*exp(ones(N_G,1)*epsi_nodes(i,:));   
            % Compute the next-period productivity levels for each integration       
            % node using the process for productivity; n_nodes-by-N
                  
            k2  =  Ord_Polynomial_N([k1 a1],D)*bk_D; 
            % Compute capital of period t+2 (chosen at t+1) using the
            % capital policy functions; n_nodes-by-N 

            C1 = (A*k1.^alpha.*a1 - k2+k1*(1-delta))*ones(N,1);
            % C is computed by summing up individual consumption, which in
            % turn, is found from the individual budget constraints; N_G-by-1

            c1 = C1*ones(1,N)/N;                 
            % Compute next-period consumption for N countries; n_nodes-by-N
      
            for j = 1:N
                Y(1:N_G,j) = Y(1:N_G,j)+weight_nodes(i,1)*beta*c1(1:N_G,j).^(-gam)./c(1:N_G,j).^(-gam).*(1-delta+alpha*A*k1(1:N_G,j).^(alpha-1).*a1(1:N_G,j)).*k1(1:N_G,j);
            end % N_G-by-N
        end
 
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
    bk_hat_D = X_G(1:N_G,:)\Y(1:N_G,:); 
                              % Compute new coefficients of the capital 
                              % policy functions using the chosen 
                              % approximation method
    bk_D = kdamp*bk_hat_D + (1-kdamp)*bk_D; 
                              % Update the coefficients of the capital  
                              % policy functions using damping 
                                     
    % 15.2.6 Store the capital series 
    %--------------------------------
    k_old = k1;       % The stored capital series will be used for checking 
                      % the convergence on the subsequent iteration

    end;

     % 15.2.7 The EDS and CGA output for the polynomial solution of degree D
     % ---------------------------------------------------------------------
     BK(1:npol(D),1:N,D) = bk_D; % Store the coefficients of the polynomial  
                                 % of degree D that approximates capital 
                                 % policy functions of N countries 
                                      
     time_EDS_CGA(D) = toc;      % Time needed to compute the polynomial  
                                 % solution of degree D 
end                         

% 16. Accuracy of the EDS_CGA solutions: errors on a stochastic simulation 
% ------------------------------------------------------------------------

% 16.1 Specify a set of points on which the accuracy is evaluated
%----------------------------------------------------------------
T_test = 10200;                    % Choose the simulation length for the  
                                   % test on a stochastic simulation, 
                                   % T_test<=10,200 

a_test = a20200(T+1:T+T_test,1:N); % Restrict the series of the productivity 
                                   % levels for the test on a stochastic 
                                   % simulation to the given T_test<=10,200  
                                   % and N<=10                             
          
k_test(1,1:N) = 1;  % Initial condition for capital (equal to steady state)

% 16.2 Choose an integration method for evaluating the accuracy of solutions
%-----------------------------------------------------------------------
IM_test = 11;                      % See paragraph 13 for the integration 
                                   % options

% To implement the test on a stochastic simulation with T_test>10,200, one
% needs to simulate new series of the productivity levels with larger T_test 
% by enabling the code in paragraph 5

% 16.3 Compute errors on a stochastic simulation for the EDS_CGA polynomial 
% solution of degrees D=1,...,D_max
% -----------------------------------------------------------------------------
for D = 1:D_max 
    
    % 16.3.1 Simulate the time series solution under the given capital-
    % policy-function coefficients, BK(:,:,D) with D=1,...,D_max 
    %------------------------------------------------------------------
    bk = BK(1:npol(D),:,D);     % The vector of coefficients of the 
                                % polynomial of degree D 
    
    for t = 1:T_test-1
        X_test = Ord_Polynomial_N([k_test(t,:) a_test(t,:)],D);
        % The basis functions of a polynomial of degree D at time t
        k_test(t+1,:) = X_test*bk;
        % Compute next-period capital using bk
    end    
    
    % 16.3.2 Errors across 10,000 points on a stochastic simulation
    % -------------------------------------------------------------
    discard = 200; % Discard the first 200 observations to remove the effect
                   % of the initial conditions 
    [Errors_mean(1,D),Errors_max(1,D), time_test(1,D)] = Accuracy_Test_N(k_test,a_test,bk,D,IM_test,alpha,gam,delta,beta,A,tau,rho,vcv,discard);

    % Errors_mean    is the unit-free average absolute approximation error  
    %                across 4N+1 optimality conditions (in log10) 
    % Errors_max     is the unit-free maximum absolute approximation error   
    %                across 4N+1 optimality conditions (in log10) 
end

% 17. Display the results for the polynomial solutions of the degrees from 
% one to D_max  
% ------------------------------------------------------------------------
disp(' '); disp('           GSSA OUTPUT:'); disp(' '); 
disp('RUNNING TIME (in seconds):'); disp('');
disp('a) for computing the solution'); 
disp('      1        2         3         4        5');disp(time_EDS_CGA);
disp('b) for implementing the accuracy test'); 
disp('      1        2         3         4        5'); disp(time_test);
disp('APPROXIMATION ERRORS (log10):'); disp(''); 
disp('a) mean error across 4N+1 optimality conditions'); 
disp('      1        2         3         4        5');disp(Errors_mean)
disp('b) max error across 4N+1 optimality conditions'); 
disp('      1        2         3         4        5');disp(Errors_max)

%save Results_N time_EDS_CGA time_test Errors_mean Errors_max kdamp IM N T BK k_test a_test IM_test alpha gam delta beta A tau rho vcv discard npol D_max T_test ;
