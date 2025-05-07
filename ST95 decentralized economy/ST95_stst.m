% Stockman and Tesar (1995), Tastes and Technology in a Two-Country Model of the Business Cycle: Explaining International Comovements, 
% American Economic Review, Vol. 85, No. 1 (Mar., 1995), pp. 168-185

% written by Katrin Rabitsch, March 2009

% ST95_stst.m

% This file finds the non-stochastic steady state using csolve, calls file
% (ST95_stst_equations.m) containing the system of equations

function [a,as,kT,kNT,kTs,kNTs,aT,aNT,aTs,aNTs,kTp,kNTp,kTsp,kNTsp,aTp,aNTp,aTsp,aNTsp,c1,c2,c1s,c2s,d,ds,cc,ccs,nT,nNT,nTs,nNTs,CC,CCs,LL,LLs,yT,yNT,yTs,yNTs,iT,iNT,iTs,iNTs,c1p,c2p,c1sp,c2sp,dp,dsp,ccp,ccsp,nTp,nNTp,nTsp,nNTsp,CCp,CCsp,LLp,LLsp,yTp,yNTp,yTsp,yNTsp,iTp,iNTp,iTsp,iNTsp,wT,wNT,wTs,wNTs,rT,rNT,rTs,rNTs,p1,p2,pT,pTs,pNT,pNTs,RER,wTp,wNTp,wTsp,wNTsp,rTp,rNTp,rTsp,rNTsp,p1p,p2p,pTp,pTsp,pNTp,pNTsp,RERp]=ST95_stst;

global i_IRorSIM

% assign parameters
[Nbar,Nsbar,aTbar,aTsbar,aNTbar,aNTsbar,omeg,sig,mu,mus,theta,thetas,alphaT,alphaNT,alphaTs,alphaNTs,delta,betta,gam,eta,rho_aTaT,rho_aTaNT,rho_aTaTs,rho_aTaNTs,rho_aNTaT,rho_aNTaNT,rho_aNTaTs,rho_aNTaNTs,rho_aTsaT,rho_aTsaNT,rho_aTsaTs,rho_aTsaNTs,rho_aNTsaT,rho_aNTsaNT,rho_aNTsaTs,rho_aNTsaNTs,eta11,eta12,eta13,eta14,eta21,eta22,eta23,eta24,eta31,eta32,eta33,eta34,eta41,eta42,eta43,eta44]=ST95_param;
% vector of starting values
x0=ones(35,1)*.1; x0(29,1)=.99; 
% calling csolve
[SS,rc]=csolve(@ST95_stst_equations,x0,[],0.00001,1000)

% storing steady state output
c1   = SS(1);          
c2   = SS(2);          
c1s  = SS(3);                
c2s  = SS(4);          
d    = SS(5);           
ds   = SS(6);          
kT   = SS(7);            
kTs  = SS(8);          
kNT  = SS(9);          
kNTs = SS(10);           
p1   = SS(11);           
p2   = SS(12);         
pNT  = SS(13);         
pNTs = SS(14);         
pT   = SS(15);         
pTs  = SS(16);           
CC   = SS(17);              
CCs  = SS(18);          
wT   = SS(19);          
wTs  = SS(20);              
wNT  = SS(21);            
wNTs = SS(22);               
rT   = SS(23);         
rTs  = SS(24);          
rNT  = SS(25);            
rNTs = SS(26);         
RER  = SS(29);         
cc   = SS(30);         
ccs  = SS(31);            
nn   = Nbar;
nns  = Nsbar;
a    = SS(34);         
as   = SS(35);        
nT   = SS(32);           
nTs  = SS(33);               
nNT  = SS(27);             
nNTs = SS(28);               
aT   = aTbar;
aTs  = aTsbar;
aNT  = aNTbar;
aNTs = aNTsbar;
iT   = ( gam*kT   - (1-delta)*kT  );
iNT  = ( gam*kNT  - (1-delta)*kNT );
iTs  = ( gam*kTs  - (1-delta)*kTs );
iNTs = ( gam*kNTs - (1-delta)*kNTs);
LL   = (1-nT -nNT);   
LLs  = (1-nTs-nNTs);  
yT   = (aT  *kT  ^alphaT  *nT  ^(1-alphaT));
yNT  = (aNT *kNT ^alphaNT *nNT ^(1-alphaNT));
yTs  = (aTs *kTs ^alphaTs *nTs ^(1-alphaTs));
yNTs = (aNTs*kNTs^alphaNTs*nNTs^(1-alphaNTs));

% displaying steady state output
if i_IRorSIM==0;
fprintf('\n\n  STEADY STATE VALUES \n');
if rc==0;     fprintf('\n  (normal solution) \n\n') 
elseif rc==4; fprintf('\n  (WARNING: maximum number of iterations reached) \n\n')
else          fprintf('\n  (WARNING: no solution) \n\n') 
end
fprintf('     c1 = %8.5f\n',c1);
fprintf('     c2 = %8.5f\n',c2);   
fprintf('    c1s = %8.5f\n',c1s);
fprintf('    c2s = %8.5f\n',c2s);  
fprintf('      d = %8.5f\n',d);  
fprintf('     ds = %8.5f\n',ds);
fprintf('     kT = %8.5f\n',kT);
fprintf('    kTs = %8.5f\n',kTs);  
fprintf('    kNT = %8.5f\n',kNT); 
fprintf('   kNTs = %8.5f\n',kNTs);
fprintf('     p1 = %8.5f\n',p1);    
fprintf('     p2 = %8.5f\n',p2);       
fprintf('    pNT = %8.5f\n',pNT);      
fprintf('   pNTs = %8.5f\n',pNTs);       
fprintf('     pT = %8.5f\n',pT);      
fprintf('    pTs = %8.5f\n',pTs);   
fprintf('     CC = %8.5f\n',CC); 
fprintf('    CCs = %8.5f\n',CCs); 
fprintf('     wT = %8.5f\n',wT);     
fprintf('    wTs = %8.5f\n',wTs);  
fprintf('    wNT = %8.5f\n',wNT);   
fprintf('   wNTs = %8.5f\n',wNTs); 
fprintf('     rT = %8.5f\n',rT);      
fprintf('    rTs = %8.5f\n',rTs);      
fprintf('    rNT = %8.5f\n',rNT);  
fprintf('   rNTs = %8.5f\n',rNTs);       
fprintf('    RER = %8.5f\n',RER);       
fprintf('     cc = %8.5f\n',cc);      
fprintf('    ccs = %8.5f\n',ccs);    
fprintf('      a = %8.5f\n',a);      
fprintf('     as = %8.5f\n',as);       
fprintf('     nT = %8.5f\n',nT);    
fprintf('    nTs = %8.5f\n',nTs); 
fprintf('    nNT = %8.5f\n',nNT);  
fprintf('   nNTs = %8.5f\n',nNTs); 
end % if i_IRorSIM
         

% Applying logs and defining 'prime'-variables
c1   = log(c1);      
c2   = log(c2);       
c1s  = log(c1s);      
c2s  = log(c2s);       
d    = log(d);      
ds   = log(ds);
kT   = log(kT);  
kTs  = log(kTs);  
kNT  = log(kNT);  
kNTs = log(kNTs);  
p1   = log(p1);      
p2   = log(p2);       
pNT  = log(pNT);      
pNTs = log(pNTs);       
pT   = log(pT);      
pTs  = log(pTs);       
CC   = log(CC);      
CCs  = log(CCs);       
wT   = log(wT);      
wTs  = log(wTs);       
wNT  = log(wNT);      
wNTs = log(wNTs);       
rT   = log(rT);      
rTs  = log(rTs);       
rNT  = log(rNT);      
rNTs = log(rNTs);       
RER  = log(RER);       
cc   = log(cc);      
ccs  = log(ccs);       
nn   = log(Nbar);
nns  = log(Nsbar);
nT   = log(nT);      
nTs  = log(nTs);       
nNT  = log(nNT);      
nNTs = log(nNTs);  
aT   = log(aTbar);
aTs  = log(aTsbar);
aNT  = log(aNTbar);
aNTs = log(aNTsbar);
iT   = log(iT);
iNT  = log(iNT);
iTs  = log(iTs);
iNTs = log(iNTs);
LL   = log(LL);
LLs  = log(LLs);
yT   = log(yT);
yNT  = log(yNT);
yTs  = log(yTs);
yNTs = log(yNTs);


c1p   = (c1);      
c2p   = (c2);       
c1sp  = (c1s);      
c2sp  = (c2s);       
dp    = (d);      
dsp   = (ds);
kTp   = (kT);  
kTsp  = (kTs);  
kNTp  = (kNT);  
kNTsp = (kNTs);  
p1p   = (p1);      
p2p   = (p2);       
pNTp  = (pNT);      
pNTsp = (pNTs);       
pTp   = (pT);      
pTsp  = (pTs);       
CCp   = (CC);      
CCsp  = (CCs);       
wTp   = (wT);      
wTsp  = (wTs);       
wNTp  = (wNT);      
wNTsp = (wNTs);       
rTp   = (rT);      
rTsp  = (rTs);       
rNTp  = (rNT);      
rNTsp = (rNTs);       
RERp  = (RER);       
ccp   = (cc);      
ccsp  = (ccs);       
nnp   = (Nbar);
nnsp  = (Nsbar);
nTp   = (nT);      
nTsp  = (nTs);       
nNTp  = (nNT);      
nNTsp = (nNTs);  
aTp   = (aT);
aTsp  = (aTs);
aNTp  = (aNT);
aNTsp = (aNTs);
iTp   = (iT);
iNTp  = (iNT);
iTsp  = (iTs);
iNTsp = (iNTs);
LLp   = (LL);
LLsp  = (LLs);
yTp   = (yT);
yNTp  = (yNT);
yTsp  = (yTs);
yNTsp = (yNTs);







%%
% Stockman and Tesar (1995), Tastes and Technology in a Two-Country Model of the Business Cycle: Explaining International Comovements, 
% American Economic Review, Vol. 85, No. 1 (Mar., 1995), pp. 168-185

% written by Katrin Rabitsch, March 2009

% ST95_stst_equations.m

function y=ST95_stst_equations(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assign parameters
[Nbar,Nsbar,aTbar,aTsbar,aNTbar,aNTsbar,omeg,sig,mu,mus,theta,thetas,alphaT,alphaNT,alphaTs,alphaNTs,delta,betta,gam,eta,rho_aTaT,rho_aTaNT,rho_aTaTs,rho_aTaNTs,rho_aNTaT,rho_aNTaNT,rho_aNTaTs,rho_aNTaNTs,rho_aTsaT,rho_aTsaNT,rho_aTsaTs,rho_aTsaNTs,rho_aNTsaT,rho_aNTsaNT,rho_aNTsaTs,rho_aNTsaNTs,eta11,eta12,eta13,eta14,eta21,eta22,eta23,eta24,eta31,eta32,eta33,eta34,eta41,eta42,eta43,eta44]=ST95_param;
[rows,cols]=size(x);
j=1;
while j<=cols;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%variables in the system:

c1      = x(1,j);      
c2      = x(2,j);      
c1s     = x(3,j);     
c2s     = x(4,j);     
d       = x(5,j); 
ds      = x(6,j);
kT      = x(7,j);   kTp     = kT;
kTs     = x(8,j);   kTsp    = kTs;
kNT     = x(9,j);   kNTp    = kNT;
kNTs    = x(10,j);  kNTsp   = kNTs;
p1      = x(11,j);  p1p     = p1;
p2      = x(12,j);  p2p     = p2;
pNT     = x(13,j);  pNTp    = pNT;
pNTs    = x(14,j);  pNTsp   = pNTs;
pT      = x(15,j);  pTp     = pT;
pTs     = x(16,j);  pTsp    = pTs;
CC      = x(17,j);  CCp     = CC;
CCs     = x(18,j);  CCsp    = CCs;
wT      = x(19,j);  wTp     = wT;
wTs     = x(20,j);  wTsp    = wTs;
wNT     = x(21,j);  wNTp    = wNT;
wNTs    = x(22,j);  wNTsp   = wNTs;
rT      = x(23,j);  rTp     = rT;
rTs     = x(24,j);  rTsp    = rTs;
rNT     = x(25,j);  rNTp    = rNT;
rNTs    = x(26,j);  rNTsp   = rNTs;
RER     = x(29,j);  RERp    = RER;
cc      = x(30,j);      
ccs     = x(31,j);      
a       = x(34,j); 
as      = x(35,j);
nT      = x(32,j);  nTp     = nT;
nTs     = x(33,j);  nTsp    = nTs;
nNT     = x(27,j);  nNTp    = nNT;
nNTs    = x(28,j);  nNTsp   = nNTs;

% stst productivity and leisure:
aT   = 1; aTp   = aT;
aTs  = 1; aTsp  = aTs;
aNT  = 1; aNTp  = aNT;
aNTs = 1; aNTsp = aNTs;
nn  = Nbar;
nns = Nsbar;
LL = (1-nT -nNT);  LLp =LL;  
LLs= (1-nTs-nNTs); LLsp=LLs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consumption CES aggregators and demand functions:
y(1,j)  = CC -(cc ^(-mu) +d ^(-mu)) ^(-1/mu);   
y(2,j)  = CCs-(ccs^(-mus)+ds^(-mus))^(-1/mus);   

y(3,j)  = cc -(pT  ^(-1/(1+mu)) * CC);   
y(4,j)  = ccs-(pTs ^(-1/(1+mus))* CCs);  
y(5,j)  = d  -(pNT ^(-1/(1+mu)) * CC);   
y(6,j)  = ds- (pNTs^(-1/(1+mus))* CCs);  

y(7,j)  = c1 -(   theta  *(p1/pT)       ^(-1) * cc);   
y(8,j)  = c1s-(   thetas *(p1/(pTs*RER))^(-1) * ccs);  
y(9,j)  = c2 -((1-theta) *(p2/pT)       ^(-1) * cc);   
y(10,j) = c2s-((1-thetas)*(p2/(pTs*RER))^(-1) * ccs);  
y(11,j) = cc      -(c1 ^theta *c2 ^(1-theta)) ;  
y(12,j) = ccs     -(c1s^thetas*c2s^(1-thetas));   

% Households' Intratemporal Conditions:
y(32,j) = (1/(1-sig))*CC  * a *(1-nn) ^(-1) -  wT  ;  
y(33,j) = (1/(1-sig))*CC  * a *(1-nn) ^(-1) -  wNT ;  
y(17,j) = (1/(1-sig))*CCs * as*(1-nns)^(-1) -  wTs ; 
y(18,j) = (1/(1-sig))*CCs * as*(1-nns)^(-1) -  wNTs; 

y(34,j) = nn - (nT + nNT);
y(35,j) = nn - (nT + nNT);

% Households' Intertemporal Conditions:
y(13,j) = gam*p1        - (betta*(p1p        *(1-delta)+ rTp  )); 
y(14,j) = gam*(p2/RER)  - (betta*((p2p/RERp) *(1-delta)+ rTsp ));
y(15,j) = gam*pNT       - (betta*(pNTp       *(1-delta)+ rNTp )); 
y(16,j) = gam*pNTs      - (betta*(pNTsp      *(1-delta)+ rNTsp));  

y(19,j) = RER - (CCs^(-sig)*LLs^as)/(CC^(-sig)*LL^a);  

% Firms' Optimality Conditions and Production Functions:
y(20,j) =  wT      /p1  - (1-alphaT)  *aT  *kT  ^ alphaT  *nT  ^(-alphaT);
y(21,j) = (wTs*RER)/p2  - (1-alphaTs) *aTs *kTs ^ alphaTs *nTs ^(-alphaTs);  
y(22,j) =  wNT     /pNT - (1-alphaNT) *aNT *kNT ^ alphaNT *nNT ^(-alphaNT); 
y(23,j) =  wNTs    /pNTs- (1-alphaNTs)*aNTs*kNTs^ alphaNTs*nNTs^(-alphaNTs); 
y(24,j) =  rT      /p1  -    alphaT   *aT  *kT  ^(alphaT-1)  *nT  ^(1-alphaT);  
y(25,j) = (rTs*RER)/p2  -    alphaTs  *aTs *kTs ^(alphaTs-1) *nTs ^(1-alphaTs);  
y(26,j) =  rNT     /pNT -    alphaNT  *aNT *kNT ^(alphaNT-1) *nNT ^(1-alphaNT);  
y(27,j) =  rNTs    /pNTs-    alphaNTs *aNTs*kNTs^(alphaNTs-1)*nNTs^(1-alphaNTs); 

% Resource constraints:
y(28,j) = c1  + c1s + gam*kTp   - (1-delta)*kT   - (aT  *kT  ^alphaT  *nT  ^(1-alphaT));  
y(30,j) = c2  + c2s + gam*kTsp  - (1-delta)*kTs  - (aTs *kTs ^alphaTs *nTs ^(1-alphaTs)); 
y(29,j) =       d   + gam*kNTp  - (1-delta)*kNT  - (aNT *kNT ^alphaNT *nNT ^(1-alphaNT)); 
y(31,j) =       ds  + gam*kNTsp - (1-delta)*kNTs - (aNTs*kNTs^alphaNTs*nNTs^(1-alphaNTs));

j=j+1;
end;
