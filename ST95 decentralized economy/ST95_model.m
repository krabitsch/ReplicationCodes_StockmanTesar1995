% Stockman and Tesar (1995), Tastes and Technology in a Two-Country Model of the Business Cycle: Explaining International Comovements, 
% American Economic Review, Vol. 85, No. 1 (Mar., 1995), pp. 168-185

% written by Katrin Rabitsch, March 2009

function [nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,f,eta] = ST95_model(approx);

global i_IRorSIM

%Define parameters
syms omeg sig mu mus theta thetas a as alphaT alphaNT alphaTs alphaNTs delta betta gam eta
syms rho_aTaT   rho_aTaNT   rho_aTaTs   rho_aTaNTs   aTbar 
syms rho_aNTaT  rho_aNTaNT  rho_aNTaTs  rho_aNTaNTs  aNTbar
syms rho_aTsaT  rho_aTsaNT  rho_aTsaTs  rho_aTsaNTs  aTsbar
syms rho_aNTsaT rho_aNTsaNT rho_aNTsaTs rho_aNTsaNTs aNTsbar
syms eta11 eta12 eta13 eta14 eta21 eta22 eta23 eta24 eta31 eta32 eta33 eta34 eta41 eta42 eta43 eta44    


%Define variables 
syms kT  kNT  kTs  kNTs  aT  aNT  aTs  aNTs
syms kTp kNTp kTsp kNTsp aTp aNTp aTsp aNTsp
syms c1  c2  c1s  c2s  d  ds  cc  ccs  nT  nNT  nTs  nNTs  lamT  lamNT  lamTs  lamNTs  CC  CCs  LL  LLs  yT  yNT  yTs  yNTs  iT  iNT  iTs  iNTs 
syms c1p c2p c1sp c2sp dp dsp ccp ccsp nTp nNTp nTsp nNTsp lamTp lamNTp lamTsp lamNTsp CCp CCsp LLp LLsp yTp yNTp yTsp yNTsp iTp iNTp iTsp iNTsp
syms wT  wNT  wTs  wNTs  rT  rNT  rTs  rNTs  p1  p2  pT  pTs  pNT  pNTs  RER 
syms wTp wNTp wTsp wNTsp rTp rNTp rTsp rNTsp p1p p2p pTp pTsp pNTp pNTsp RERp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write equations, ei=1:

% Consumption CES aggregators and demand functions:
e1  = CC -(cc ^(-mu) +d ^(-mu)) ^(-1/mu); 
e2  = CCs-(ccs^(-mus)+ds^(-mus))^(-1/mus); 

e3  = cc -(pT  ^(-1/(1+mu)) * CC);  
e4  = ccs-(pTs ^(-1/(1+mus))* CCs);  
e5  = d  -(pNT ^(-1/(1+mu)) * CC);   
e6  = ds- (pNTs^(-1/(1+mus))* CCs);  

e7  = c1 -(   theta  *(p1/pT)       ^(-1) * cc);   
e8  = c1s-(   thetas *(p1/(pTs*RER))^(-1) * ccs); 
e9  = c2 -((1-theta )*(p2/pT)       ^(-1) * cc);  
e10 = c2s-((1-thetas)*(p2/(pTs*RER))^(-1) * ccs); 
e11 = cc      -(c1 ^theta *c2 ^(1-theta));  
e12 = ccs     -(c1s^thetas*c2s^(1-thetas)); 

% Households' Other Intratemporal Conditions:
e13 = (1/(1-sig))*CC  * a *LL ^(-1) -  wT  ; 
e14 = (1/(1-sig))*CC  * a *LL ^(-1) -  wNT ; 
e15 = (1/(1-sig))*CCs * as*LLs^(-1) -  wTs ; 
e16 = (1/(1-sig))*CCs * as*LLs^(-1) -  wNTs; 

% Households' Intertemporal Conditions:
e17 = gam* p1     *(CC ^(-sig)*LL ^a ) - (betta*(CCp ^(-sig)*LLp ^a )*( p1p       *(1-delta)+ rTp  )); 
e18 = gam*(p2/RER)*(CCs^(-sig)*LLs^as) - (betta*(CCsp^(-sig)*LLsp^as)*((p2p/RERp) *(1-delta)+ rTsp ));  
e19 = gam*pNT     *(CC ^(-sig)*LL ^a ) - (betta*(CCp ^(-sig)*LLp ^a )*(pNTp       *(1-delta)+ rNTp ));  
e20 = gam*pNTs    *(CCs^(-sig)*LLs^as) - (betta*(CCsp^(-sig)*LLsp^as)*(pNTsp      *(1-delta)+ rNTsp));  

% Risk Sharing equation
e21 = RER - (CCs^(-sig)*LLs^as)/(CC^(-sig)*LL^a);  

% Investment and Leisure:
e22 = iT   - ( gam*kTp   - (1-delta)*kT  );
e23 = iTs  - ( gam*kTsp  - (1-delta)*kTs );
e24 = iNT  - ( gam*kNTp  - (1-delta)*kNT );
e25 = iNTs - ( gam*kNTsp - (1-delta)*kNTs);

e26 = LL - (1-nT -nNT);  
e27 = LLs- (1-nTs-nNTs); 

% Firms' Optimality Conditions and Production Functions:
e28 =  wT      /p1  - (1-alphaT)  *aT   *kT  ^ alphaT     *nT  ^(-alphaT);
e29 = (wTs*RER)/p2  - (1-alphaTs) *aTs  *kTs ^ alphaTs    *nTs ^(-alphaTs); 
e30 =  wNT     /pNT - (1-alphaNT) *aNT  *kNT ^ alphaNT    *nNT ^(-alphaNT);  
e31 =  wNTs    /pNTs- (1-alphaNTs)*aNTs *kNTs^ alphaNTs   *nNTs^(-alphaNTs); 
e32 =  rT      /p1  -    alphaT   *aT   *kT  ^(alphaT-1)  *nT  ^(1-alphaT);   
e33 = (rTs*RER)/p2  -    alphaTs  *aTs  *kTs ^(alphaTs-1) *nTs ^(1-alphaTs);  
e34 =  rNT     /pNT -    alphaNT  *aNT  *kNT ^(alphaNT-1) *nNT ^(1-alphaNT);   
e35 =  rNTs    /pNTs-    alphaNTs *aNTs *kNTs^(alphaNTs-1)*nNTs^(1-alphaNTs); 

e36 = yT   - (aT  *kT  ^alphaT  *nT  ^(1-alphaT));
e37 = yTs  - (aTs *kTs ^alphaTs *nTs ^(1-alphaTs));
e38 = yNT  - (aNT *kNT ^alphaNT *nNT ^(1-alphaNT));
e39 = yNTs - (aNTs*kNTs^alphaNTs*nNTs^(1-alphaNTs));

% Resource constraints:
e40 = c1  + c1s + gam*kTp   - (1-delta)*kT   - (aT  *kT  ^alphaT  *nT  ^(1-alphaT));  
e41 = c2  + c2s + gam*kTsp  - (1-delta)*kTs  - (aTs *kTs ^alphaTs *nTs ^(1-alphaTs)); 
e42 =       d   + gam*kNTp  - (1-delta)*kNT  - (aNT *kNT ^alphaNT *nNT ^(1-alphaNT));
e43 =       ds  + gam*kNTsp - (1-delta)*kNTs - (aNTs*kNTs^alphaNTs*nNTs^(1-alphaNTs));

% Exogenous processes:
e44 = aTp  -(rho_aTaT  *aT +rho_aTaNT  *aNT +rho_aTaTs  *aTs +rho_aTaNTs  *aNTs +(1-rho_aTaT  - rho_aTaNT  - rho_aTaTs  - rho_aTaNTs)  *aTbar);  
e45 = aNTp -(rho_aNTaT *aT +rho_aNTaNT *aNT +rho_aNTaTs *aTs +rho_aNTaNTs *aNTs +(1-rho_aNTaT - rho_aNTaNT - rho_aNTaTs - rho_aNTaNTs) *aNTbar);  
e46 = aTsp -(rho_aTsaT *aT +rho_aTsaNT *aNT +rho_aTsaTs *aTs +rho_aTsaNTs *aNTs +(1-rho_aTsaT - rho_aTsaNT - rho_aTsaTs - rho_aTsaNTs) *aTsbar);  
e47 = aNTsp-(rho_aNTsaT*aT +rho_aNTsaNT*aNT +rho_aNTsaTs*aTs +rho_aNTsaNTs*aNTs +(1-rho_aNTsaT- rho_aNTsaNT- rho_aNTsaTs- rho_aNTsaNTs)*aNTsbar);  

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

f = [e1;e2;e3;e4;e5;e6;e7;e8;e9;e10;e11;e12;e13;e14;e15;e16;e17;e18;e19;e20;e21;e22;e23;e24;e25;e26;e27;e28;e29;e30;e31;e32;e33;e34;e35;e36;e37;e38;e39;e40;e41;e42;e43;e44;e45;e46;e47];

% Define the vector of controls, y, and states, x
x  = [kT  kTs  kNT  kNTs  aT  aTs  aNT  aNTs];
xp = [kTp kTsp kNTp kNTsp aTp aTsp aNTp aNTsp];
y  = [c1  c2  c1s  c2s  d  ds  nT  nNT  nTs  nNTs  yT  yNT  yTs  yNTs  iT  iNT  iTs  iNTs...
      cc  ccs  CC  CCs  LL  LLs  wT  wNT  wTs  wNTs  rT  rNT  rTs  rNTs  p1  p2  pT  pTs  pNT  pNTs  RER];  
yp = [c1p c2p c1sp c2sp dp dsp nTp nNTp nTsp nNTsp yTp yNTp yTsp yNTsp iTp iNTp iTsp iNTsp...
      ccp ccsp CCp CCsp LLp LLsp wTp wNTp wTsp wNTsp rTp rNTp rTsp rNTsp p1p p2p pTp pTsp pNTp pNTsp RERp]; 

%Make f a function of the logarithm of the state and control vector
f = subs(f, [x,y,xp,yp], exp([x,y,xp,yp]));

%Compute analytical derivatives of f
[nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx]=anal_deriv(f,x,y,xp,yp,approx);
