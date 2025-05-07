% Stockman and Tesar (1995), Tastes and Technology in a Two-Country Model
% of the Business Cycle: Explaining International Comovements, American
% Economic Review, Vol. 85, No. 1 (Mar., 1995), pp. 168-185

% written by Katrin Rabitsch, March 2009

function [nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,f,eta] = ST95_planner_model(approx);

global i_IRorSIM

%Define parameters
syms omeg sig mu mus theta thetas a as alphaT alphaNT alphaTs alphaNTs delta betta gam eta
syms rho_aTaT rho_aTaNT rho_aTaTs rho_aTaNTs aTbar 
syms rho_aNTaT rho_aNTaNT rho_aNTaTs rho_aNTaNTs aNTbar
syms rho_aTsaT rho_aTsaNT rho_aTsaTs rho_aTsaNTs aTsbar
syms rho_aNTsaT rho_aNTsaNT rho_aNTsaTs rho_aNTsaNTs aNTsbar
syms eta11 eta12 eta13 eta14 eta21 eta22 eta23 eta24 eta31 eta32 eta33 eta34 eta41 eta42 eta43 eta44    

%Define variables 
syms kT  kNT  kTs  kNTs  aT  aNT  aTs  aNTs
syms kTp kNTp kTsp kNTsp aTp aNTp aTsp aNTsp
syms c1  c2  c1s  c2s  d  ds  nT  nNT  nTs  nNTs  lamT  lamNT  lamTs  lamNTs  CC  CCs  LL  LLs  yT  yNT  yTs  yNTs  iT  iNT  iTs  iNTs  RER  cc  ccs
syms c1p c2p c1sp c2sp dp dsp nTp nNTp nTsp nNTsp lamTp lamNTp lamTsp lamNTsp CCp CCsp LLp LLsp yTp yNTp yTsp yNTsp iTp iNTp iTsp iNTsp RERp ccp ccsp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % derive first order conditions
% syms L
% 
% L  =    omeg * (1/(1-sig))*(( (c1 ^theta *c2 ^(1-theta)) ^(-mu) + d ^(-mu)) ^(1/(-mu)) )^(1-sig)*(1-nT -nNT )^a  + ...
%      (1-omeg)* (1/(1-sig))*(( (c1s^thetas*c2s^(1-thetas))^(-mus)+ ds^(-mus))^(1/(-mus)))^(1-sig)*(1-nTs-nNTs)^as - ...
%     lamT  *(c1  + c2  + gam*kTp   - (1-delta)*kT   - (aT  *kT  ^alphaT  *nT  ^(1-alphaT)))  - ...
%     lamTs *(c1s + c2s + gam*kTsp  - (1-delta)*kTs  - (aTs *kTs ^alphaTs *nTs ^(1-alphaTs))) - ...
%     lamNT *(d   + gam*kNTp  - (1-delta)*kNT  - (aNT *kNT ^alphaNT *nNT ^(1-alphaNT)))       - ...
%     lamNTs*(ds  + gam*kNTsp - (1-delta)*kNTs - (aNTs*kNTs^alphaNTs*nNTs^(1-alphaNTs)))
% 
% diff(L,'c1')
% diff(L,'c2')
% diff(L,'c1s')
% diff(L,'c2s')
% diff(L,'d')
% diff(L,'ds')
% diff(L,'nT')
% diff(L,'nNT')
% diff(L,'nTs')
% diff(L,'nNTs')

% FOCs wrt tradable consumption, c1, c2, c1s, c2s:
e1  =    omeg * CC ^(1-sig+mu)  *LL ^a  *(c1 ^theta *c2 ^(1-theta)) ^(-mu)  *     theta  *c1^(-1)  - lamT;  
e2  =    omeg * CC ^(1-sig+mu)  *LL ^a  *(c1 ^theta *c2 ^(1-theta)) ^(-mu)  *  (1-theta) *c2^(-1)  - lamTs; 
e3  = (1-omeg)* CCs^(1-sig+mus) *LLs^as *(c1s^thetas*c2s^(1-thetas))^(-mus) *    thetas  *c1s^(-1) - lamT;  
e4  = (1-omeg)* CCs^(1-sig+mus) *LLs^as *(c1s^thetas*c2s^(1-thetas))^(-mus) * (1-thetas) *c2s^(-1) - lamTs;

% FOCs wrt nontradable consumption, d, ds:
e5  =    omeg * CC ^(1-sig+mu)  *LL ^a  * d ^(-mu-1)  - lamNT; 
e6  = (1-omeg)* CCs^(1-sig+mus) *LLs^as * ds^(-mus-1) - lamNTs; 

% definitions of tradable consumption baskets
e36 = cc      -(c1 ^theta *c2 ^(1-theta));  
e37 = ccs     -(c1s^thetas*c2s^(1-thetas)); 

% FOCs wrt labor:
e7  =    omeg * (1/(1-sig))*CC ^(1-sig) * a *LL ^(a-1)  -  lamT  * (1-alphaT)  *aT  *kT  ^alphaT  *nT  ^(-alphaT);  
e8  = (1-omeg)* (1/(1-sig))*CCs^(1-sig) * as*LLs^(as-1) -  lamTs * (1-alphaTs) *aTs *kTs ^alphaTs *nTs ^(-alphaTs); 
e9  =    omeg * (1/(1-sig))*CC ^(1-sig) * a *LL ^(a-1)  -  lamNT * (1-alphaNT) *aNT *kNT ^alphaNT *nNT ^(-alphaNT); 
e10 = (1-omeg)* (1/(1-sig))*CCs^(1-sig) * as*LLs^(as-1) -  lamNTs* (1-alphaNTs)*aNTs*kNTs^alphaNTs*nNTs^(-alphaNTs);

% FOCs wrt capital:
e11 = lamT  *gam - (betta*lamTp  *(1-delta + alphaT  *aTp  *kTp  ^(alphaT-1)  *nTp  ^(1-alphaT)));  
e13 = lamTs *gam - (betta*lamTsp *(1-delta + alphaTs *aTsp *kTsp ^(alphaTs-1) *nTsp ^(1-alphaTs))); 
e12 = lamNT *gam - (betta*lamNTp *(1-delta + alphaNT *aNTp *kNTp ^(alphaNT-1) *nNTp ^(1-alphaNT))); 
e14 = lamNTs*gam - (betta*lamNTsp*(1-delta + alphaNTs*aNTsp*kNTsp^(alphaNTs-1)*nNTsp^(1-alphaNTs)));

% Resource constraints:
e15 = c1 + c1s  + gam*kTp   - (1-delta)*kT   - (aT  *kT  ^alphaT  *nT  ^(1-alphaT));  
e16 = c2 + c2s  + gam*kTsp  - (1-delta)*kTs  - (aTs *kTs ^alphaTs *nTs ^(1-alphaTs)); 
e17 =       d   + gam*kNTp  - (1-delta)*kNT  - (aNT *kNT ^alphaNT *nNT ^(1-alphaNT));
e18 =       ds  + gam*kNTsp - (1-delta)*kNTs - (aNTs*kNTs^alphaNTs*nNTs^(1-alphaNTs)); 

% Consumption indices and leisure:
e19 = CC  - ( (c1 ^theta *c2 ^(1-theta)) ^(-mu) + d ^(-mu)) ^(1/(-mu));  
e20 = CCs - ( (c1s^thetas*c2s^(1-thetas))^(-mus)+ ds^(-mus))^(1/(-mus)); 
e21 = LL - (1-nT -nNT);   
e22 = LLs- (1-nTs-nNTs); 

% Output and investment:
e27 = yT   - (aT  *kT  ^alphaT  *nT  ^(1-alphaT));
e28 = yNT  - (aNT *kNT ^alphaNT *nNT ^(1-alphaNT));
e29 = yTs  - (aTs *kTs ^alphaTs *nTs ^(1-alphaTs));
e30 = yNTs - (aNTs*kNTs^alphaNTs*nNTs^(1-alphaNTs));
e31 = iT   - ( gam*kTp   - (1-delta)*kT  );
e32 = iNT  - ( gam*kNTp  - (1-delta)*kNT );
e33 = iTs  - ( gam*kTsp  - (1-delta)*kTs );
e34 = iNTs - ( gam*kNTsp - (1-delta)*kNTs);

e35 = RER - (CCs^(-sig)*LLs^as)/(CC^(-sig)*LL^a);

% Exogenous processes:
e23 = aTp  -(rho_aTaT  *aT +rho_aTaNT  *aNT +rho_aTaTs  *aTs +rho_aTaNTs  *aNTs +(1-rho_aTaT  - rho_aTaNT  - rho_aTaTs  - rho_aTaNTs)  *aTbar);  
e24 = aNTp -(rho_aNTaT *aT +rho_aNTaNT *aNT +rho_aNTaTs *aTs +rho_aNTaNTs *aNTs +(1-rho_aNTaT - rho_aNTaNT - rho_aNTaTs - rho_aNTaNTs) *aNTbar);  
e25 = aTsp -(rho_aTsaT *aT +rho_aTsaNT *aNT +rho_aTsaTs *aTs +rho_aTsaNTs *aNTs +(1-rho_aTsaT - rho_aTsaNT - rho_aTsaTs - rho_aTsaNTs) *aTsbar);  
e26 = aNTsp-(rho_aNTsaT*aT +rho_aNTsaNT*aNT +rho_aNTsaTs*aTs +rho_aNTsaNTs*aNTs +(1-rho_aNTsaT- rho_aNTsaNT- rho_aNTsaTs- rho_aNTsaNTs)*aNTsbar);  

eta = [  0     0     0     0  ;
         0     0     0     0  ;
         0     0     0     0  ;
         0     0     0     0  ;
       eta11 eta12 eta13 eta14;
       eta21 eta22 eta23 eta24;
       eta31 eta32 eta33 eta34;
       eta41 eta42 eta43 eta44];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE FUNCTION f, DEFINE CONTROLS AND STATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = [e1;e2;e3;e4;e5;e6;e7;e8;e9;e10;e11;e12;e13;e14;e15;e16;e17;e18;e19;e20;e21;e22;e23;e24;e25;e26;e27;e28;e29;e30;e31;e32;e33;e34;e35;e36;e37];

% Define the vector of controls, y, and states, x
x  = [kT  kTs  kNT  kNTs  aT  aTs  aNT  aNTs];
xp = [kTp kTsp kNTp kNTsp aTp aTsp aNTp aNTsp];
y  = [c1  c2  c1s  c2s  d  ds  nT  nNT  nTs  nNTs  yT  yNT  yTs  yNTs  iT  iNT  iTs  iNTs  cc  ccs  CC  CCs  LL  LLs  lamT  lamNT  lamTs  lamNTs  RER];
yp = [c1p c2p c1sp c2sp dp dsp nTp nNTp nTsp nNTsp yTp yNTp yTsp yNTsp iTp iNTp iTsp iNTsp ccp ccsp CCp CCsp LLp LLsp lamTp lamNTp lamTsp lamNTsp RERp];

%Make f a function of the logarithm of the state and control vector
f = subs(f, [x,y,xp,yp], exp([x,y,xp,yp]));

%Compute analytical derivatives of f
[nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx]=anal_deriv(f,x,y,xp,yp,approx);
