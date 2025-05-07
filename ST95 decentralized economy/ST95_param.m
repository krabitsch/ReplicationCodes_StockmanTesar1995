% Stockman and Tesar (1995), Tastes and Technology in a Two-Country Model of the Business Cycle: Explaining International Comovements, 
% American Economic Review, Vol. 85, No. 1 (Mar., 1995), pp. 168-185

% written by Katrin Rabitsch, March 2009

% ST95_param.m

% This file sets the parameter values

function [Nbar,Nsbar,aTbar,aTsbar,aNTbar,aNTsbar,omeg,sig,mu,mus,theta,thetas,alphaT,alphaNT,alphaTs,alphaNTs,delta,betta,gam,eta,rho_aTaT,rho_aTaNT,rho_aTaTs,rho_aTaNTs,rho_aNTaT,rho_aNTaNT,rho_aNTaTs,rho_aNTaNTs,rho_aTsaT,rho_aTsaNT,rho_aTsaTs,rho_aTsaNTs,rho_aNTsaT,rho_aNTsaNT,rho_aNTsaTs,rho_aNTsaNTs,eta11,eta12,eta13,eta14,eta21,eta22,eta23,eta24,eta31,eta32,eta33,eta34,eta41,eta42,eta43,eta44]=ST95_param;

global i_IRorSIM


omeg        = 0.5;           % weight on Home country in social planner's problem
sig         = 2;             % intertemporal elasticity of substitution
elast_mu    = .44;           % elasticity of substitution between tradables and nontradables
mu          = 1/elast_mu-1;  
mus         = mu;       
theta       = 0.5;           % Home's weight on domestic tradables in tradables-basket
thetas      = 1-theta;       % Foreign's weight on domestic tradables in tradables-basket
% a           = 1/(-0.315);  % parameter that makes households work 20% (Nbar) of their time endowment at steady state
% as          = a;
alphaT      = 1-0.61;        % capital share in tradables-production, Home
alphaNT     = 1-0.56;        % capital share in non-tradables-production, Home
alphaTs     = alphaT;        % capital share in tradables-production, Foreign
alphaNTs    = alphaNT;       % capital share in non-tradables-production, Foreign
delta       = .1;            % rate of depreciation of capital stock
betta       = 0.96;          % discount factor
gam         = 1.0273;        % annual growth rate

% Shock persistences for the vector of shocks [aT, aNT, aTs, aNTs]
if     i_IRorSIM==0; % shock persistences used for Impulse Responses
    OMEG= [ 0.9   0     0    0;
            0     0.9   0    0;
            0     0     0.9  0;
            0     0     0    0.9];    
elseif i_IRorSIM==1; % shock persistences used for simulation
    OMEG= [ 0.154  0.040  -0.199  0.262;
           -0.150  0.632  -0.110  0.125;
           -0.199  0.262   0.154  0.040;
           -0.110  0.125  -0.015  0.632];
end 

rho_aTaT  =OMEG(1,1); rho_aTaNT  =OMEG(1,2); rho_aTaTs  =OMEG(1,3); rho_aTaNTs  =OMEG(1,4);
rho_aNTaT =OMEG(2,1); rho_aNTaNT =OMEG(2,2); rho_aNTaTs =OMEG(2,3); rho_aNTaNTs =OMEG(2,4);
rho_aTsaT =OMEG(3,1); rho_aTsaNT =OMEG(3,2); rho_aTsaTs =OMEG(3,3); rho_aTsaNTs =OMEG(3,4);
rho_aNTsaT=OMEG(4,1); rho_aNTsaNT=OMEG(4,2); rho_aNTsaTs=OMEG(4,3); rho_aNTsaNTs=OMEG(4,4);

% Variance-Covariance matrix and definition of eta
 VCV=[ 3.62 1.23 1.21 0.51;
       1.23 1.99 0.51 0.27;
       1.21 0.51 3.62 1.23;       
       0.51 0.27 1.23 1.99];
eta         = [zeros(4,4); chol(VCV)/100];

eta11=eta(1,1); eta12=eta(1,2); eta13=eta(1,3); eta14=eta(1,4); 
eta21=eta(2,1); eta22=eta(2,2); eta23=eta(2,3); eta24=eta(2,4); 
eta31=eta(3,1); eta32=eta(3,2); eta33=eta(3,3); eta34=eta(3,4); 
eta41=eta(4,1); eta42=eta(4,2); eta43=eta(4,3); eta44=eta(4,4); 

% steady state parameters
Nbar=1/5;                    % Home steady state labor (as ration of total time endowment)
Nsbar=1/5;                   % Foreign steady state labor (as ration of total time endowment)
aTbar=1;                     % steady state value of productivity
aTsbar=1;                    % steady state value of productivity
aNTbar=1;                    % steady state value of productivity 
aNTsbar=1;                   % steady state value of productivity