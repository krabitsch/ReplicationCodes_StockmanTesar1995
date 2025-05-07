% Stockman and Tesar (1995), Tastes and Technology in a Two-Country Model of the Business Cycle: Explaining International Comovements, 
% American Economic Review, Vol. 85, No. 1 (Mar., 1995), pp. 168-185

% written by Katrin Rabitsch, March 2009

% ST95_planner_stst.m

% This file finds the non-stochastic steady state using csolve, calls file
% (ST95_planner_stst_equations.m) containing the system of equations

function [a,as,c1,c2,c1s,c2s,d,ds,kT,kTs,kNT,kNTs,lamT,lamTs,lamNT,lamNTs,aT,aTs,aNT,aNTs,nT,nTs,nNT,nNTs,CC,CCs,LL,LLs,c1p,c2p,c1sp,c2sp,dp,dsp,kTp,kTsp,kNTp,kNTsp,lamTp,lamTsp,lamNTp,lamNTsp,aTp,aTsp,aNTp,aNTsp,nTp,nTsp,nNTp,nNTsp,CCp,CCsp,LLp,LLsp,iT,iNT,iTs,iNTs,RER,iTp,iNTp,iTsp,iNTsp,RERp,yT,yNT,yTs,yNTs,yTp,yNTp,yTsp,yNTsp,cc,ccs,ccp,ccsp]=ST95_planner_stst;

global i_IRorSIM
         
% assign parameters
[Nbar,Nsbar,aTbar,aTsbar,aNTbar,aNTsbar,omeg,sig,mu,mus,theta,thetas,alphaT,alphaNT,alphaTs,alphaNTs,delta,betta,gam,eta,rho_aTaT,rho_aTaNT,rho_aTaTs,rho_aTaNTs,rho_aNTaT,rho_aNTaNT,rho_aNTaTs,rho_aNTaNTs,rho_aTsaT,rho_aTsaNT,rho_aTsaTs,rho_aTsaNTs,rho_aNTsaT,rho_aNTsaNT,rho_aNTsaTs,rho_aNTsaNTs,eta11,eta12,eta13,eta14,eta21,eta22,eta23,eta24,eta31,eta32,eta33,eta34,eta41,eta42,eta43,eta44]=ST95_planner_param;
% vector of starting values
x0=ones(16,1)*.05; 
x0=[0.0747;
    0.0747;
    0.0747;
    0.0747;
    0.1074;
    0.1074;
    0.4838;
    0.4838;
    0.4140;
    0.4140;
    0.1241;
    0.1241;
    0.0759;
    0.0759;
   -.5;
   -.5];

% calling csolve
[SS,rc]=csolve(@ST95_planner_stst_equations,x0,[],0.00001,1000)

% storing steady state output
c1     = SS(1);          
c2     = SS(2);              
c1s    = SS(3);             
c2s    = SS(4);               
d      = SS(5);           
ds     = SS(6);          
kT     = SS(7);          
kTs    = SS(8);          
kNT    = SS(9);          
kNTs   = SS(10);         
nn     = Nbar;
nns    = Nsbar;
nT     = SS(11);              
nTs    = SS(12);            
nNT    = SS(13);         
nNTs   = SS(14);         
% lamT   = SS(15);         
% lamTs  = SS(16);         
% lamNT  = SS(17);           
% lamNTs = SS(18);         
a      = SS(15);         
as     = SS(16);         

CC   = ( (c1 ^theta *c2 ^(1-theta)) ^(-mu) + d ^(-mu)) ^(1/(-mu));  
CCs  = ( (c1s^thetas*c2s^(1-thetas))^(-mus)+ ds^(-mus))^(1/(-mus)); 
cc   = (c1 ^theta *c2 ^(1-theta));   
ccs  = (c1s^thetas*c2s^(1-thetas));  
LL   = (1-nT -nNT);   
LLs  = (1-nTs-nNTs);  
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
RER  = (CCs^(-sig)*LLs^as)/(CC^(-sig)*LL^a);
lamT    = (1-omeg)* CCs^(1-sig+mus) *LLs^as *(c1s^thetas*c2s^(1-thetas))^(-mus) *     thetas *c1s^(-1) ; 
lamTs   = (1-omeg)* CCs^(1-sig+mus) *LLs^as *(c1s^thetas*c2s^(1-thetas))^(-mus) * (1-thetas) *c2s^(-1);
lamNT   = omeg    * CC ^(1-sig+mu)  *LL ^a  * d ^(-mu-1)  ;   
lamNTs  = (1-omeg)* CCs^(1-sig+mus) *LLs^as * ds^(-mus-1); 


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
fprintf('     nT = %8.5f\n',nT); 
fprintf('    nTs = %8.5f\n',nTs);    
fprintf('    nNT = %8.5f\n',nNT);      
fprintf('   nNTs = %8.5f\n',nNTs);       
fprintf('   lamT = %8.5f\n',lamT);      
fprintf('  lamTs = %8.5f\n',lamTs);       
fprintf('  lamNT = %8.5f\n',lamNT);    
fprintf(' lamNTs = %8.5f\n',lamNTs);       
fprintf('      a = %8.5f\n',a);      
fprintf('     as = %8.5f\n',as);       
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
CC   = log(CC);      
CCs  = log(CCs);       
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
lamT   = log(lamT);
lamNT  = log(lamNT);
lamTs  = log(lamTs);
lamNTs = log(lamNTs);
RER    = log(RER);


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
CCp   = (CC);      
CCsp  = (CCs);       
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
lamTp   = (lamT);
lamNTp  = (lamNT);
lamTsp  = (lamTs);
lamNTsp = (lamNTs);
RERp    = (RER);







%%
% Stockman and Tesar (1995), Tastes and Technology in a Two-Country Model of the Business Cycle: Explaining International Comovements, 
% American Economic Review, Vol. 85, No. 1 (Mar., 1995), pp. 168-185

% written by Katrin Rabitsch, March 2009

% ST95_planner_stst_equations.m

function y=ST95_planner_stst_equations(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assign parameters
[Nbar,Nsbar,aTbar,aTsbar,aNTbar,aNTsbar,omeg,sig,mu,mus,theta,thetas,alphaT,alphaNT,alphaTs,alphaNTs,delta,betta,gam,eta,rho_aTaT,rho_aTaNT,rho_aTaTs,rho_aTaNTs,rho_aNTaT,rho_aNTaNT,rho_aNTaTs,rho_aNTaNTs,rho_aTsaT,rho_aTsaNT,rho_aTsaTs,rho_aTsaNTs,rho_aNTsaT,rho_aNTsaNT,rho_aNTsaTs,rho_aNTsaNTs,eta11,eta12,eta13,eta14,eta21,eta22,eta23,eta24,eta31,eta32,eta33,eta34,eta41,eta42,eta43,eta44]=ST95_planner_param;
[rows,cols]=size(x);
j=1;
while j<=cols;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% % variables in the system:
c1      = x(1,j);      
c2      = x(2,j);      
c1s     = x(3,j);     
c2s     = x(4,j);     
d       = x(5,j); 
ds      = x(6,j);
kT      = x(7,j);   kTp     = x(7,j);
kTs     = x(8,j);   kTsp    = x(8,j);
kNT     = x(9,j);   kNTp    = x(9,j);
kNTs    = x(10,j);  kNTsp   = x(10,j);
nT      = x(11,j);  nTp     = x(11,j);
nTs     = x(12,j);  nTsp    = x(12,j);
nNT     = x(13,j);  nNTp    = x(13,j);
nNTs    = x(14,j);  nNTsp   = x(14,j);
% lamT    = x(15,j);  lamTp   = x(15,j);
% lamTs   = x(16,j);  lamTsp  = x(16,j);
% lamNT   = x(17,j);  lamNTp  = x(17,j);
% lamNTs  = x(18,j);  lamNTsp = x(18,j);
a       = x(15,j); 
as      = x(16,j);

% stst productivity:
aT   = 1; aTp   = aT;
aTs  = 1; aTsp  = aTs;
aNT  = 1; aNTp  = aNT;
aNTs = 1; aNTsp = aNTs;
nn =Nbar;
nns=Nsbar;

% Consumption indices and leisure:
CC  = ( (c1 ^theta *c2 ^(1-theta)) ^(-mu) + d ^(-mu)) ^(1/(-mu));  
CCs = ( (c1s^thetas*c2s^(1-thetas))^(-mus)+ ds^(-mus))^(1/(-mus)); 
LL  = (1-nT -nNT);   
LLs = (1-nTs-nNTs);  

lamT    = (1-omeg)* CCs^(1-sig+mus) *LLs^as *(c1s^thetas*c2s^(1-thetas))^(-mus) *     thetas *c1s^(-1) ; 
lamTs   = (1-omeg)* CCs^(1-sig+mus) *LLs^as *(c1s^thetas*c2s^(1-thetas))^(-mus) * (1-thetas) *c2s^(-1);
lamNT   = omeg    * CC ^(1-sig+mu)  *LL ^a  * d ^(-mu-1)  ;   
lamNTs  = (1-omeg)* CCs^(1-sig+mus) *LLs^as * ds^(-mus-1); 
lamTp   = lamT;
lamTsp  = lamTs;
lamNTp  = lamNT;
lamNTsp = lamNTs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM:
% FOCs wrt tradable consumption, c1, c2, c1s, c2s:
y(1,j)  = omeg    * CC ^(1-sig+mu)  *LL ^a  *(c1 ^theta *c2 ^(1-theta)) ^(-mu)  *     theta  *c1 ^(-1) - lamT;  
y(2,j)  = omeg    * CC ^(1-sig+mu)  *LL ^a  *(c1 ^theta *c2 ^(1-theta)) ^(-mu)  * (1-theta)  *c2 ^(-1) - lamTs;
% y(3,j)  = (1-omeg)* CCs^(1-sig+mus) *LLs^as *(c1s^thetas*c2s^(1-thetas))^(-mus) *     thetas *c1s^(-1) - lamT; 
% y(4,j)  = (1-omeg)* CCs^(1-sig+mus) *LLs^as *(c1s^thetas*c2s^(1-thetas))^(-mus) * (1-thetas) *c2s^(-1) - lamTs;

% FOCs wrt nontradable consumption, d, ds:
% y(5,j)  = omeg    * CC ^(1-sig+mu)  *LL ^a  * d ^(-mu-1)  - lamNT;   
% y(6,j)  = (1-omeg)* CCs^(1-sig+mus) *LLs^as * ds^(-mus-1) - lamNTs; 

% % FOCs wrt labor:
y(7,j)  =    omeg * (1/(1-sig))*CC ^(1-sig) * a *LL ^(a-1)  -  lamT  * (1-alphaT)  *aT  *kT  ^alphaT  *nT  ^(-alphaT);  
y(8,j)  = (1-omeg)* (1/(1-sig))*CCs^(1-sig) * as*LLs^(as-1) -  lamTs * (1-alphaTs) *aTs *kTs ^alphaTs *nTs ^(-alphaTs); 
y(9,j)  =    omeg * (1/(1-sig))*CC ^(1-sig) * a *LL ^(a-1)  -  lamNT * (1-alphaNT) *aNT *kNT ^alphaNT *nNT ^(-alphaNT); 
y(10,j) = (1-omeg)* (1/(1-sig))*CCs^(1-sig) * as*LLs^(as-1) -  lamNTs* (1-alphaNTs)*aNTs*kNTs^alphaNTs*nNTs^(-alphaNTs); 

% FOCs wrt capital:
y(11,j) = lamT  *gam - (betta*lamTp  *(1-delta + alphaT  *aTp  *kTp  ^(alphaT-1)  *nTp  ^(1-alphaT)));  
y(13,j) = lamTs *gam - (betta*lamTsp *(1-delta + alphaTs *aTsp *kTsp ^(alphaTs-1) *nTsp ^(1-alphaTs))); 
y(12,j) = lamNT *gam - (betta*lamNTp *(1-delta + alphaNT *aNTp *kNTp ^(alphaNT-1) *nNTp ^(1-alphaNT))); 
y(14,j) = lamNTs*gam - (betta*lamNTsp*(1-delta + alphaNTs*aNTsp*kNTsp^(alphaNTs-1)*nNTsp^(1-alphaNTs)));

% Resource constraints:
y(3,j)  = c1 + c1s + gam*kTp   - (1-delta)*kT   - (aT  *kT  ^alphaT  *nT  ^(1-alphaT));  
y(4,j)  = c2 + c2s + gam*kTsp  - (1-delta)*kTs  - (aTs *kTs ^alphaTs *nTs ^(1-alphaTs)); 
y(5,j)  =       d  + gam*kNTp  - (1-delta)*kNT  - (aNT *kNT ^alphaNT *nNT ^(1-alphaNT));  
y(6,j)  =      ds  + gam*kNTsp - (1-delta)*kNTs - (aNTs*kNTs^alphaNTs*nNTs^(1-alphaNTs));

y(15,j) = nn - (nT + nNT);
y(16,j) = nn - (nT + nNT);

j=j+1;
end;
