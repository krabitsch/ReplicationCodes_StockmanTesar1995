% Stockman and Tesar (1995), Tastes and Technology in a Two-Country Model of the Business Cycle: Explaining International Comovements, 
% American Economic Review, Vol. 85, No. 1 (Mar., 1995), pp. 168-185

% written by Katrin Rabitsch, March 2009

% Main file; calls the model file, gets solution, plots impulse responses and simulates the model

clear
global i_IRorSIM

approx    = 1
i_IRorSIM = 0;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compute numerical derivatives of model
[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f,eta]=ST95_model(approx); 
if approx==1;
anal_deriv_print2f('ST95',fx,fxp,fy,fyp,f,eta);   
elseif approx==2;
anal_deriv_print2f('ST95',fx,fxp,fy,fyp,f,eta, fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx);
end
%Numerical evaluation
%Assign values to parameters and steady-state varia1bles
[Nbar,Nsbar,aTbar,aTsbar,aNTbar,aNTsbar,omeg,sig,mu,mus,theta,thetas,alphaT,alphaNT,alphaTs,alphaNTs,delta,betta,gam,eta,rho_aTaT,rho_aTaNT,rho_aTaTs,rho_aTaNTs,rho_aNTaT,rho_aNTaNT,rho_aNTaTs,rho_aNTaNTs,rho_aTsaT,rho_aTsaNT,rho_aTsaTs,rho_aTsaNTs,rho_aNTsaT,rho_aNTsaNT,rho_aNTsaTs,rho_aNTsaNTs,eta11,eta12,eta13,eta14,eta21,eta22,eta23,eta24,eta31,eta32,eta33,eta34,eta41,eta42,eta43,eta44]=ST95_param;
[a,as,kT,kNT,kTs,kNTs,aT,aNT,aTs,aNTs,kTp,kNTp,kTsp,kNTsp,aTp,aNTp,aTsp,aNTsp,c1,c2,c1s,c2s,d,ds,cc,ccs,nT,nNT,nTs,nNTs,CC,CCs,LL,LLs,yT,yNT,yTs,yNTs,iT,iNT,iTs,iNTs,c1p,c2p,c1sp,c2sp,dp,dsp,ccp,ccsp,nTp,nNTp,nTsp,nNTsp,CCp,CCsp,LLp,LLsp,yTp,yNTp,yTsp,yNTsp,iTp,iNTp,iTsp,iNTsp,wT,wNT,wTs,wNTs,rT,rNT,rTs,rNTs,p1,p2,pT,pTs,pNT,pNTs,RER,wTp,wNTp,wTsp,wNTsp,rTp,rNTp,rTsp,rNTsp,p1p,p2p,pTp,pTsp,pNTp,pNTsp,RERp]=ST95_stst;
ST95_num_eval;
% Get solution to first order approximation
[gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp)
% Get solution to second order approximation
if approx==2;
%Second-order approximation
[gxx,hxx] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx); 
[gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta);    
end % if approx






%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impulse Responses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ic1=1;      ic2=2;     ic1s=3;     ic2s=4;      id=5;       ids=6;     inT=7;      inNT=8;    inTs=9;    inNTs=10;
iyT=11;     iyNT=12;   iyTs=13;    iyNTs=14;    iiT=15;     iiNT=16;   iiTs=17;    iiNTs=18;  icc=19;    iccs=20;
iCC=21;     iCCs=22;   iLL=23;     iLLs=24;     iwT=25;     iwNT=26;   iwTs=27;    iwNTs=28;  irT=29;    irNT=30;   
irTs=31;    irNTs=32;  ip1=33;     ip2=34;      ipT=35;     ipTs=36;   ipNT=37;    ipNTs=38;  iRER=39;
% lagged variables
ikT=size(gx,1)+1;   ikTs=size(gx,1)+2;  ikNT=size(gx,1)+3;  ikNTs=size(gx,1)+4;

T=50;           % number of IR periods                
kk=[1:T];

% shock to domestic Nontraded-Good Productivity
x0=zeros(size(hx,1),1);  x0(end-1,1)=1;    
IR=ir(gx,hx,x0,T); 
figure % Figure 4a) and 4b) in Cleveland Federal Reserve Working Paper version
plot(kk,IR(1:T,ic1) ,'b +',kk,IR(1:T,ic2) ,'b :',kk,IR(1:T,id) ,'b *',kk,IR(1:T,inT) ,'b -',kk,IR(1:T,inNT) ,'b --',kk,IR(1:T,ikT) ,'b -.',kk,IR(1:T,ikNT) ,'b o') 
legend('c_1','c_2','d','n_T','n_{NT}','k_T','k_{NT}'),
title('Home-Country Response to Nontraded-Good Productivity Shock')
figure
plot(kk,IR(1:T,ic1s),'b +',kk,IR(1:T,ic2s),'b :',kk,IR(1:T,ids),'b *',kk,IR(1:T,inTs),'b -',kk,IR(1:T,inNTs),'b --',kk,IR(1:T,ikTs),'b -.',kk,IR(1:T,ikNTs),'b o') 
legend('c_1^{\ast}','c_2^{\ast}','d^{\ast}','n_T^{\ast}','n_{NT}^{\ast}','k_T^{\ast}','k_{NT}^{\ast}'),
title('Foreign-Country Response to Nontraded-Good Productivity Shock')

% shock to domestic Traded-Good Productivity
x0=zeros(size(hx,1),1);  x0(end-3,1)=1;    
IR=ir(gx,hx,x0,T); 
kk=[1:T];
figure % Figure 3a) and 3b) in Cleveland Federal Reserve Working Paper version
plot(kk,IR(1:T,ic1) ,'b +',kk,IR(1:T,ic2) ,'b :',kk,IR(1:T,id) ,'b *',kk,IR(1:T,inT) ,'b -',kk,IR(1:T,inNT) ,'b --',kk,IR(1:T,ikT) ,'b -.',kk,IR(1:T,ikNT) ,'b o') 
legend('c_1','c_2','d','n_T','n_{NT}','k_T','k_{NT}'),
title('Home-Country Response to Traded-Good Productivity Shock')
figure
plot(kk,IR(1:T,ic1s),'b +',kk,IR(1:T,ic2s),'b :',kk,IR(1:T,ids),'b *',kk,IR(1:T,inTs),'b -',kk,IR(1:T,inNTs),'b --',kk,IR(1:T,ikTs),'b -.',kk,IR(1:T,ikNTs),'b o') 
legend('c_1^{\ast}','c_2^{\ast}','d^{\ast}','n_T^{\ast}','n_{NT}^{\ast}','k_T^{\ast}','k_{NT}^{\ast}'),
title('Foreign-Country Response to Traded-Good Productivity Shock')







%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if approx==2;
i_IRorSIM=1;
[Nbar,Nsbar,aTbar,aTsbar,aNTbar,aNTsbar,omeg,sig,mu,mus,theta,thetas,alphaT,alphaNT,alphaTs,alphaNTs,delta,betta,gam,eta,rho_aTaT,rho_aTaNT,rho_aTaTs,rho_aTaNTs,rho_aNTaT,rho_aNTaNT,rho_aNTaTs,rho_aNTaNTs,rho_aTsaT,rho_aTsaNT,rho_aTsaTs,rho_aTsaNTs,rho_aNTsaT,rho_aNTsaNT,rho_aNTsaTs,rho_aNTsaNTs]=ST95_param;
[a,as,kT,kNT,kTs,kNTs,aT,aNT,aTs,aNTs,kTp,kNTp,kTsp,kNTsp,aTp,aNTp,aTsp,aNTsp,c1,c2,c1s,c2s,d,ds,cc,ccs,nT,nNT,nTs,nNTs,CC,CCs,LL,LLs,yT,yNT,yTs,yNTs,iT,iNT,iTs,iNTs,c1p,c2p,c1sp,c2sp,dp,dsp,ccp,ccsp,nTp,nNTp,nTsp,nNTsp,CCp,CCsp,LLp,LLsp,yTp,yNTp,yTsp,yNTsp,iTp,iNTp,iTsp,iNTsp,wT,wNT,wTs,wNTs,rT,rNT,rTs,rNTs,p1,p2,pT,pTs,pNT,pNTs,RER,wTp,wNTp,wTsp,wNTsp,rTp,rNTp,rTsp,rNTsp,p1p,p2p,pTp,pTsp,pNTp,pNTsp,RERp]=ST95_stst;
num_eval;
% Get solution to first order approximation
[gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);
% Get solution to second order approximation
if approx==2;
%Second-order approximation
[gxx,hxx] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx); 
[gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta);    
end % if approx

simT=5000
SIG=1;
% e=zeros(simT,4);
e=randn(simT,4);
x0=zeros(8,1)';
[Y,X] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, SIG, x0, e);
%  Y ==> controls and definitions
%  X ==> states

    %compute variances and covariances
     %(throw away first 50 values)
     Z=[Y X];
     Z=Z(51:end,:);
     %applying HP Filter
     % lambda=100;
     % [Zhp]=hpfilter(Z,lambda);
     % or on unfiltered series
     Zhp=Z;
     %get standard deviation in percent
     varcovmatrix=cov(Zhp);
     stdevs=((diag(varcovmatrix))'.^(1/2))*100;

disp('Some Moments of Stockman and Tesar (AER 1995), Table 6: Technology shocks only')     
disp('Standard deviations, traded-good sector')
fprintf(1,'Output         = %10.6f\n',stdevs(1,iyT))
fprintf(1,'Labor          = %10.6f\n',stdevs(1,inT))
fprintf(1,'Investment     = %10.6f\n',stdevs(1,iiT))
fprintf(1,'Consumption    = %10.6f\n',stdevs(1,icc))

disp('Standard deviations, nontraded-good sector')
fprintf(1,'Output         = %10.6f\n',stdevs(1,iyNT))
fprintf(1,'Labor          = %10.6f\n',stdevs(1,inNT))
fprintf(1,'Investment     = %10.6f\n',stdevs(1,iiNT))
fprintf(1,'Consumption    = %10.6f\n',stdevs(1,id))
     
end % if approx==2