<<<<<<< HEAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Main Program for
%
%   STOCHASTIC SIMULATIONS OF NF-kB PATHWAY  
%   S. Tay et al. 2010 Nature
%   
%   Calls: Model, Parameters, AllCellPlotting and AvarageCellPlotting 
%
%   After running MainFile you can run 
%   AllCellPlotting and AvarageCellPlotting
%   
%   Saves all data in 'last' - can be used to make plots latter on
%
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear;              %reset all
% clc;                %clear comand window
starttime=clock;    %current time
rng(sum(1000*clock), 'twister');

%#####################################
%######### Simulation setup ##########
%#####################################
TNF_List = linspace(-3,1,20);
TNF=TNF_List(1);                       % TNF dose 

ANa=2; AN=2; ANR=2;           % ANa=2 - # IKBa alleles, AN=2 - # A20 alleles,  ANR=2 - # Reporter gene alleles
%ANa=2; AN=2; ANR=2;           % ANa=2 - # IKBa alleles, AN=2 - # A20 alleles,  ANR=2 - # Reporter gene alleles

%Set AN=0 to study A20 knockout%

N=1000;                          % number of cell to be simulated
RS=1; %RampSpeed unit/hr
fold=1.5; %Fold Change
%####################################################################
%###### Simulation time points                                #######
%### Various time protocols can be studied within this frame  #######
%####################################################################

t000=10*60*60;         %10h randomization of initial conditions
t00=10*60*60;          %10h equilibrium waiting time

t0=.1*60*60;             % 1 step, time when TNF is being introduced into the system (in seconds)

tw1=.1*60*60;             % 2 step, length of TNF stimulation
te1=2*60*60;             % 3 step length of first break   (White breaks: 3600s, 6000s, 12000s, our break 170*60s)
%  
% tw2 = RS*60*60;         % 4 step length of second TNF stimulation
% te2 = 0*60*60;         % 5 step length of second break
% 
% tw3 = RS*60*60;          % 6 step length of third TNF stimulation  
% te3 = 0*60*60;          % 7 step length of third break
% 
% tw4 = RS*60*60;          % 8 step length of 4th TNF stimulation  
% te4 = 0*60*60;          % 9 step length of 4th break
% 
% tw5 = RS*60*60;          % 10 step length of 5th TNF stimulation
% te5 = 0*60*60;          % 11 step length of 5th break
% 
% tw6 = RS*60*60;          % 12 step length of 6th TNF stimulation
% te6 = 0*60*60;          % 13 step length of 6th break
% 
% tw7 = RS*60*60;          % 14 step length of 7th TNF stimulation
% te7 = 0*60*60;          % 15 step length of 7th break
% 
% tw8 = RS*60*60;          % 16 step length of 8th TNF stimulation
% te8 = 0*60*60;          % 17 step length of 8th break
% 
% tw9 = RS*60*60;          % 18 step length of 9th TNF stimulation
% te9 = 0*60*60;          % 19 step length of 9th break
% 
% tw10 = RS*60*60;          % 20 step length of 9th TNF stimulation
% te10 = 0*60*60;          % 21 step length of 9th break

% ############################################################

tt=100;            % forword time for ODEs solving 
                     
YYY=0;                         %matrix of average, all variables y0(i)(t)
totalNFkB=0;                        %total nuclear NF-kB 
GGa=0; GG=0;  GGT=0; GGR=0;    %status of Ikba,  A20, TNF and reporter genes
Bb=0;                          %number of active receptors
MM=0;
NFF=0;

for i=1:N           %beginning the mean loop 

    %#################################
    %###### Initial conditions #######
    %#################################

    i                  %cell nummber
    
    
    [meanNFkB,par1NFkB,par2NFkB,meanTNFR1,par1TNFR1,par2TNFR1,k4,ka20,AB,kv,q1,q2,c1,c3,c4,c5,k1,k2,k3,a1,a2,a3,c1a,c5a,c6a,i1,i1a,e1a,e2a,dt,tp,KN,KNN,ka,ki,actRateTNFR1,inactRateTNFR1,Tdeg,q1r,q2r,q2rr,c1r,c1rr,c3r]=Parameters;
    
    %%%% Randomizations of total TNF receptors and NF-kB levels %%%%%%%
    
      
     NF=round(meanNFkB*exp(par2NFkB+randn*par1NFkB))                 %NF-kB level
     while NF > 10*meanNFkB
     NF=round(meanNFkB*exp(par2NFkB+randn*par1NFkB))
     end
     
%     NF=meanNFkB;     %uncomment to remove extrinsic noise 
%    
%    Lognormal distribution with Median=meanNFkB, Mean=meanNFkB*Exp(par1NFkB^2/2), %
%    Variance = meanNFkB^2 * (Exp (par1NFkB^2 -1) * Exp(par1NFkB^2)
   

     M=round(meanTNFR1*exp(par2TNFR1+randn*par1TNFR1))            % number of TNFR1 receptors
     while M > 10*meanTNFR1
     M=round(meanTNFR1*exp(par2TNFR1+randn*par1TNFR1))
     end
     
%    M=meanTNFR1;         %uncomment to remove extrinsic noise  
%    
%    Lognormal distribution with Median=meanTNFR1, Mean=meanTNFR1*Exp(par1TNFR1^2/2), %
%    Variance = meanTNFR1^2 * (Exp (par1TNFR1^2 -1) * Exp(par1TNFR1^2)
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    

    y0=zeros(1,19);     %initial conditions set to zero and next: 

    y0(14)=NF;          %NF-kB is given in cytoplasmic complex(IkBa|NFkB) (total NF-kB kept constant), standard = 10^5
    y0(2)=2*10^5;       %initial IKKn, total IKK kept constant  
    y0(11)=0.14*y0(14); %free cytoplasmic IkBa protein 
    y0(12)=0.06*y0(14); %free nuclear IkBa protein
    y0(13)=10;          %IkBa mRNA
    y0(10)=10;          %10 A20 mRNA
    y0(9)=10000;        %10000 A20 protein

    
    y0(10)=AB*y0(10);                   
    y0(9)=AB*y0(9); 

    Ga=0;               % initial status of IkBa promoter
    G=0;                % initial status of A20 promoter
    GT=0;               % initial status of TNF promoter
    GR=0;               % initial status of reporter gene promoter
    B=0;                % initial number of active receptors
    yy0=y0;             % initial conditions y0(i)   

    %###################################################################
    %###### -1 step - randomization of initial condition  #######
    %###################################################################

    realtime=0;                     %simulated time   
    phase=round(rand*t000/dt)*dt;   %random initial time (dt -simulation time step -10s)  
    tspan=[0:dt:tt];                %time for which the solution is derived to find the switching time, tt=1h
   
    while (realtime<phase)
        [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
        Yact=Y0(:,8);               %amount of NF-kBn  
        Yin=Y0(:,12);               %amount of IkBan 
        TR=Y0(:,16);                %TNF level  
        Gax=Ga;Gx=G;GRx=GR;Bx=B;
        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M) ;  %function determining the change of gene status, calls statuschange
        
        tc=T0(mk);                  %time when the status changes 
        yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration  
        realtime=realtime+tc;
    end

    if (realtime>phase)
        nn=(realtime-phase)/dt;
        yy0=Y0(mk-nn,:);
        Ga=Gax;G=Gx;GR=GRx;B=Bx;
    end                            %status before the last change it occured outside of the time interval

    clear Yact Yin Y0 T0 nn mk phase tc;

    %####################################################
    %###### 0 step - waiting for "equilibrium"    #######
    %####################################################

    realtime=0;
    
    while (realtime<t00)
        [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
        Yact=Y0(:,8);               %amount of NF-kBn  
        Yin=Y0(:,12);               %amount of IkBan  
        TR=Y0(:,16);                %TNF level  
        Gax=Ga;Gx=G;GRx=GR;Bx=B;
        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, calls statuschange
        
        tc=T0(mk);                  %time when the status changes 
        yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration  
        realtime=realtime+tc;
    end

    if (realtime>t00)
        nn=(realtime-t00)/dt;
        yy0=Y0(mk-nn,:);
        Ga=Gax;G=Gx;GR=GRx;B=Bx;
    end                            %status before the last change it occured ouside of the time interval

    clear Yact Yin Y0 T0 nn mk tc; 

    %#######################################
    %###### 1 step- still no TNF #######
    %#######################################
    
    realtime=0;
    ga=[Ga];g=[G];gR=[GR];                  %saves activity of IkBa A20 reporter genes 
    bb=[B];                                 %saves number of active receptors
    Y=yy0;                                  %variables where single cell run is stored
    T=zeros(1,1);
    
    while (realtime<t0)
        [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
        Yact=Y0(:,8);               %amount of NF-kBn  
        Yin=Y0(:,12);               %amount of IkBan  
        TR=Y0(:,16);                %TNF level  
        Gax=Ga;Gx=G;GRx=GR;Bx=B;
        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
        tc=T0(mk);                  %time when the status changes 
        yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
        Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
        T=[T;T0(2:mk)+realtime];     
        ga=[ga;Gax*ones(mk-1,1)];
        g=[g;Gx*ones(mk-1,1)];
        gR=[gR;GRx*ones(mk-1,1)];
        bb=[bb;Bx*ones(mk-1,1)];
        realtime=realtime+tc;
    end

    
        nn=(realtime-t0)/dt;
        x=size(Y);
        Y=Y(1:(x(1)-nn),:);
        T=T(1:(x(1)-nn));
        
        Y0(mk-nn,16)=0;        %setting TNF ON for the next step
        
        yy0=Y0(mk-nn,:);
        x1=length(ga);
        ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
        Ga=Gax;G=Gx;GR=GRx;B=Bx;
 
        
        tnfcaselist = [];
        for counter = 1:size(TNF_List,2)
            yy0(16)=10^TNF_List(counter);        %setting TNF ON for the next step
            tnfcaselist(counter,:,:) = TNFloop(t0,tw1,te1,tspan,yy0,Ga,G,GR,B,Y0,GT,AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M,T0,mk,Y,T,ga,g,gR,bb,dt);
        end        
        

%     %####################################
%     %######  2 step TNF on 1 time        #######
%     %####################################
%     
%   
%     realtime=t0;
% 
%     while (realtime<t0+tw1)
%        [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GTx=GT;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end
%     
%     
%         nn=(realtime-t0-tw1)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         Y0(mk-nn,16)=0;        %setting TNF OFF for the next step
%      
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GT=GTx;GR=GRx;B=Bx;
%  
%      
%     %########################################
%     %###### 3 step TNF washed out 1 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1;
%     
%      while (realtime<t0+tw1+te1)
%         [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%         [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%      end
%     
%    
%         nn=(realtime-t0-tw1-te1)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         %Y0(mk-nn,16)=TNF;        %setting TNF ON for the next step
%         Y0(mk-nn,16)=0;        %HT20170508-setting TNF ON for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
 
%%%%%%%%%%%%%%%%5 here
% 
%     %########################################
%     %###### 4 step TNF on for 2 time ###############
%     %########################################
%     
%     
%     realtime=t0+tw1+te1;
%     
%     while (realtime<t0+tw1+te1+tw2)
%         [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         Y0(mk-nn,16)=0;        %setting TNF OFF for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%   
%     %########################################
%     %###### 5 step TNF washed out 2 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2;
%     
%       while (realtime<t0+tw1+te1+tw2+te2)
%         [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%         [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         %Y0(mk-nn,16)=TNF;        %setting TNF ON for the next step
%         Y0(mk-nn,16)=TNF*fold*fold;        %HT20170508-setting TNF ON for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
% 
%   
%     %########################################
%     %###### 6 step TNF on for the 3 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2;
%     
%      while (realtime<t0+tw1+te1+tw2+te2+tw3)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         Y0(mk-nn,16)=0;        %setting TNF OFF for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%   
%     %########################################
%     %###### 7 step TNF washed out 3 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3;
%     
%       while (realtime<t0+tw1+te1+tw2+te2+tw3+te3)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         %Y0(mk-nn,16)=TNF;        %HT:setting TNF ON for the next step
%         Y0(mk-nn,16)=TNF*fold*fold*fold;        %setting TNF ON for the next step
% 
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%  
%      %########################################
%     %###### 8 step TNF on for the 4 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3;
%     
%      while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         Y0(mk-nn,16)=0;        %setting TNF OFF for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%   
%     %########################################
%     %###### 9 step TNF washed out 4 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4;
%     
%       while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         %Y0(mk-nn,16)=TNF;        %setting TNF ON for the next step
%         Y0(mk-nn,16)=TNF*fold*fold*fold*fold;        %setting TNF ON for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%         
%      %########################################
%     %###### 10 step TNF on for the 5 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4;
%     
%      while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         Y0(mk-nn,16)=0;        %setting TNF OFF for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%   
%     %########################################
%     %###### 11 step TNF washed out 5 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5;
%     
%       while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         %Y0(mk-nn,16)=TNF;        %setting TNF ON for the next step
%         Y0(mk-nn,16)=TNF*fold*fold*fold*fold*fold;        %setting TNF ON for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%         
%              %########################################
%     %###### 12 step TNF on for the 6 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5;
%     
%      while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5-tw6)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         Y0(mk-nn,16)=0;        %setting TNF OFF for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%   
%     %########################################
%     %###### 13 step TNF washed out 6 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6;
%     
%       while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5-tw6-te6)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         %Y0(mk-nn,16)=TNF;        %setting TNF ON for the next step
%         Y0(mk-nn,16)=TNF*fold*fold*fold*fold*fold*fold;        %setting TNF ON for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%         
%              %########################################
%     %###### 14 step TNF on for the 7 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6;
%     
%      while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5-tw6-te6-tw7)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         Y0(mk-nn,16)=0;        %setting TNF OFF for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%   
%     %########################################
%     %###### 15 step TNF washed out 7 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7;
%     
%       while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5-tw6-te6-tw7-te7)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         %Y0(mk-nn,16)=TNF;        %setting TNF ON for the next step
%         Y0(mk-nn,16)=TNF*fold*fold*fold*fold*fold*fold*fold;        %setting TNF ON for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%         
%              %########################################
%     %###### 16 step TNF on for the 8 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7;
%     
%      while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5-tw6-te6-tw7-te7-tw8)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         Y0(mk-nn,16)=0;        %setting TNF OFF for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%   
%     %########################################
%     %###### 17 step TNF washed out 8 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8;
%     
%       while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8+te8)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5-tw6-te6-tw7-te7-tw8-te8)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         %Y0(mk-nn,16)=TNF;        %setting TNF ON for the next step
%         Y0(mk-nn,16)=TNF*fold*fold*fold*fold*fold*fold*fold*fold;        %setting TNF ON for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%         
%              %########################################
%     %###### 18 step TNF on for the 9 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8+te8;
%     
%      while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8+te8+tw9)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5-tw6-te6-tw7-te7-tw8-te8-tw9)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         Y0(mk-nn,16)=0;        %setting TNF OFF for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%   
%     %########################################
%     %###### 19 step TNF washed out 9 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8+te8+tw9;
%     
%       while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8+te8+tw9+te9)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5-tw6-te6-tw7-te7-tw8-te8-tw9-te9)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         %Y0(mk-nn,16)=TNF;        %setting TNF ON for the next step
%         Y0(mk-nn,16)=TNF*fold*fold*fold*fold*fold*fold*fold*fold*fold;        %setting TNF ON for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%         
%         
%                      %########################################
%     %###### 20 step TNF on for the 10 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8+te8+tw9+te9;
%     
%      while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8+te8+tw9+te9+tw10)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5-tw6-te6-tw7-te7-tw8-te8-tw9-te9-tw10)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         Y0(mk-nn,16)=0;        %setting TNF OFF for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%   
%     %########################################
%     %###### 21 step TNF washed out 9 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8+te8+tw9+te9+tw10;
%     
%       while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8+te8+tw9+te9+te9+tw10+te10)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5-tw6-te6-tw7-te7-tw8-te8-tw9-te9-tw10-te10)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         %Y0(mk-nn,16)=TNF;        %setting TNF ON for the next step
%         Y0(mk-nn,16)=0;        %setting TNF ON for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%         
%         
%         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%        
       

%     YYY=YYY+Y;
%     GGa=GGa+ga;
%     GG=GG+g;
%     GGR=GGR+gR;
%     Bb=Bb+bb;
    
MM(i)=M;
NFF(i)=NF;
    
%#############################
%###### DATA FOR PLOTS #######
%#############################

    
    XXX(i,:,:,:)=Y; 
%     XB(i,:)=bb(:);   
%     XG(i,:)=g(:); 
%     XGa(i,:)=ga(:);    
%     XGR(i,:)=gR(:);  % data for all cells 
    
end                % end of the Main LOOP

% clear tspan tspan1 g ga T0 Y0 Yact Yin Y y0 yy0 tindex;
% 
% 
% 
% YYY=YYY/N;                       %average over population
% GGa=GGa/N;
% GG=GG/N;
% GGR=GGR/N;
% Bb=Bb/N;
% T=(T-t0)/60;

 
% % Data Processing : Counts Cells Responding to First and Second Pulse
% % designated for two pulses only
% 
% RecetorsAverage=sum(MM)/N
% NFaverage=sum(NFF)/N
% 
% NFKB=XXX(:,:,8)+XXX(:,:,15);
% both=0;
% onlyfirst=0;
% onlysecond=0;
% 
% both2=0;
% onlyfirst2=0;
% onlysecond2=0;
% 
% for i=1:N
%     s=0;
%     f=0;
%     s2=0;
%     f2=0;
%     aa=NFKB(i,:);
%     c1=size(aa);
%     c=round(c1(2)/2);
%     a1=aa(i:c);
%     a2=aa(c:2*c-1);
%     
% if max(a1)>meanNFkB/10  
%     f=1;
% end    
% if max(a1)>meanNFkB/5  
%     f2=1;
% end    
% if max(a2)>meanNFkB/10
%     s=1;
% end  
% if max(a2)>meanNFkB/5
%     s2=1;
% end 
%     both=both+f*s;
%     onlyfirst=onlyfirst+f-f*s;
%     onlysecond=onlysecond+s-f*s;
%     
%     both2=both2+f2*s2;
%     onlyfirst2=onlyfirst2+f2-f2*s2;
%     onlysecond2=onlysecond2+s2-f2*s2;
% end
% 
% Threshold=0.1
% 
% both
% onlyfirst
% onlysecond
% any=both+onlyfirst+onlysecond
% 
% Threshold=0.2  % Used for this study
% 
% both2
% onlyfirst2
% onlysecond2
% any2=both2+onlyfirst2+onlysecond2
% 
% %% End of data processing %%%


simulation_time=etime(clock,starttime) %simulation time seconds
save last

% %HT 2017_0508
% mcol=sqrt(N);
% mrow=sqrt(N);
% ind=1:1:N;
% X_limm=RS*60*11;
% yrightlim=TNF*fold^9;
% 
% figure(1);hold on
% for mm=1:length(ind)
% %     subplot(mrow, mcol, mm);
%     plot(T, XXX(mm,:,8))
% %     set(gca, 'xlim', [0,1000], 'ylim', [0,10*10^4], 'yscale', 'linear');
%     title('NFkB response', 'fontsize', 8);
% end
%  figure(2);clf;
% for mm=1:length(ind)
% %     subplot(mrow, mcol, mm);
%     plot(T, XXX(mm,:,16))
%     set(gca, 'xlim', [0,1000], 'ylim', [0,20], 'yscale', 'linear');
%     title('TNF conc.', 'fontsize', 8);
% end
%  
% figure(3);clf;
% for mm=1:length(ind)
% %     subplot(mrow, mcol, mm);
%     yyaxis left
%     plot(T, XXX(mm,:,8),'b')
%     set(gca, 'xlim', [0,X_limm], 'yscale', 'linear');
%     %set(gca, 'xlim', [0,X_limm], 'ylim', [0,10*10^4], 'yscale', 'linear');
%     ylabel('Nuc. P65')
%     hold on
%     yyaxis right
%     plot(T, XXX(mm,:,16),'r')
%     set(gca, 'xlim', [0,X_limm], 'ylim', [0,yrightlim], 'yscale', 'linear');
%     ylabel('Ligand Conc.')
% end
%AllCellPlotting
=======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Main Program for
%
%   STOCHASTIC SIMULATIONS OF NF-kB PATHWAY  
%   S. Tay et al. 2010 Nature
%   
%   Calls: Model, Parameters, AllCellPlotting and AvarageCellPlotting 
%
%   After running MainFile you can run 
%   AllCellPlotting and AvarageCellPlotting
%   
%   Saves all data in 'last' - can be used to make plots latter on
%
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear;              %reset all
% clc;                %clear comand window
starttime=clock;    %current time
rng(sum(1000*clock), 'twister');

%#####################################
%######### Simulation setup ##########
%#####################################
TNF_List = linspace(-3,1,20);
TNF=TNF_List(1);                       % TNF dose 

ANa=2; AN=2; ANR=2;           % ANa=2 - # IKBa alleles, AN=2 - # A20 alleles,  ANR=2 - # Reporter gene alleles
%ANa=2; AN=2; ANR=2;           % ANa=2 - # IKBa alleles, AN=2 - # A20 alleles,  ANR=2 - # Reporter gene alleles

%Set AN=0 to study A20 knockout%

N=1000;                          % number of cell to be simulated
RS=1; %RampSpeed unit/hr
fold=1.5; %Fold Change
%####################################################################
%###### Simulation time points                                #######
%### Various time protocols can be studied within this frame  #######
%####################################################################

t000=10*60*60;         %10h randomization of initial conditions
t00=10*60*60;          %10h equilibrium waiting time

t0=.1*60*60;             % 1 step, time when TNF is being introduced into the system (in seconds)

tw1=.1*60*60;             % 2 step, length of TNF stimulation
te1=2*60*60;             % 3 step length of first break   (White breaks: 3600s, 6000s, 12000s, our break 170*60s)
%  
% tw2 = RS*60*60;         % 4 step length of second TNF stimulation
% te2 = 0*60*60;         % 5 step length of second break
% 
% tw3 = RS*60*60;          % 6 step length of third TNF stimulation  
% te3 = 0*60*60;          % 7 step length of third break
% 
% tw4 = RS*60*60;          % 8 step length of 4th TNF stimulation  
% te4 = 0*60*60;          % 9 step length of 4th break
% 
% tw5 = RS*60*60;          % 10 step length of 5th TNF stimulation
% te5 = 0*60*60;          % 11 step length of 5th break
% 
% tw6 = RS*60*60;          % 12 step length of 6th TNF stimulation
% te6 = 0*60*60;          % 13 step length of 6th break
% 
% tw7 = RS*60*60;          % 14 step length of 7th TNF stimulation
% te7 = 0*60*60;          % 15 step length of 7th break
% 
% tw8 = RS*60*60;          % 16 step length of 8th TNF stimulation
% te8 = 0*60*60;          % 17 step length of 8th break
% 
% tw9 = RS*60*60;          % 18 step length of 9th TNF stimulation
% te9 = 0*60*60;          % 19 step length of 9th break
% 
% tw10 = RS*60*60;          % 20 step length of 9th TNF stimulation
% te10 = 0*60*60;          % 21 step length of 9th break

% ############################################################

tt=100;            % forword time for ODEs solving 
                     
YYY=0;                         %matrix of average, all variables y0(i)(t)
totalNFkB=0;                        %total nuclear NF-kB 
GGa=0; GG=0;  GGT=0; GGR=0;    %status of Ikba,  A20, TNF and reporter genes
Bb=0;                          %number of active receptors
MM=0;
NFF=0;

for i=1:N           %beginning the mean loop 

    %#################################
    %###### Initial conditions #######
    %#################################

    i                  %cell nummber
    
    
    [meanNFkB,par1NFkB,par2NFkB,meanTNFR1,par1TNFR1,par2TNFR1,k4,ka20,AB,kv,q1,q2,c1,c3,c4,c5,k1,k2,k3,a1,a2,a3,c1a,c5a,c6a,i1,i1a,e1a,e2a,dt,tp,KN,KNN,ka,ki,actRateTNFR1,inactRateTNFR1,Tdeg,q1r,q2r,q2rr,c1r,c1rr,c3r]=Parameters;
    
    %%%% Randomizations of total TNF receptors and NF-kB levels %%%%%%%
    
      
     NF=round(meanNFkB*exp(par2NFkB+randn*par1NFkB))                 %NF-kB level
     while NF > 10*meanNFkB
     NF=round(meanNFkB*exp(par2NFkB+randn*par1NFkB))
     end
     
%     NF=meanNFkB;     %uncomment to remove extrinsic noise 
%    
%    Lognormal distribution with Median=meanNFkB, Mean=meanNFkB*Exp(par1NFkB^2/2), %
%    Variance = meanNFkB^2 * (Exp (par1NFkB^2 -1) * Exp(par1NFkB^2)
   

     M=round(meanTNFR1*exp(par2TNFR1+randn*par1TNFR1))            % number of TNFR1 receptors
     while M > 10*meanTNFR1
     M=round(meanTNFR1*exp(par2TNFR1+randn*par1TNFR1))
     end
     
%    M=meanTNFR1;         %uncomment to remove extrinsic noise  
%    
%    Lognormal distribution with Median=meanTNFR1, Mean=meanTNFR1*Exp(par1TNFR1^2/2), %
%    Variance = meanTNFR1^2 * (Exp (par1TNFR1^2 -1) * Exp(par1TNFR1^2)
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    

    y0=zeros(1,19);     %initial conditions set to zero and next: 

    y0(14)=NF;          %NF-kB is given in cytoplasmic complex(IkBa|NFkB) (total NF-kB kept constant), standard = 10^5
    y0(2)=2*10^5;       %initial IKKn, total IKK kept constant  
    y0(11)=0.14*y0(14); %free cytoplasmic IkBa protein 
    y0(12)=0.06*y0(14); %free nuclear IkBa protein
    y0(13)=10;          %IkBa mRNA
    y0(10)=10;          %10 A20 mRNA
    y0(9)=10000;        %10000 A20 protein

    
    y0(10)=AB*y0(10);                   
    y0(9)=AB*y0(9); 

    Ga=0;               % initial status of IkBa promoter
    G=0;                % initial status of A20 promoter
    GT=0;               % initial status of TNF promoter
    GR=0;               % initial status of reporter gene promoter
    B=0;                % initial number of active receptors
    yy0=y0;             % initial conditions y0(i)   

    %###################################################################
    %###### -1 step - randomization of initial condition  #######
    %###################################################################

    realtime=0;                     %simulated time   
    phase=round(rand*t000/dt)*dt;   %random initial time (dt -simulation time step -10s)  
    tspan=[0:dt:tt];                %time for which the solution is derived to find the switching time, tt=1h
   
    while (realtime<phase)
        [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
        Yact=Y0(:,8);               %amount of NF-kBn  
        Yin=Y0(:,12);               %amount of IkBan 
        TR=Y0(:,16);                %TNF level  
        Gax=Ga;Gx=G;GRx=GR;Bx=B;
        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M) ;  %function determining the change of gene status, calls statuschange
        
        tc=T0(mk);                  %time when the status changes 
        yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration  
        realtime=realtime+tc;
    end

    if (realtime>phase)
        nn=(realtime-phase)/dt;
        yy0=Y0(mk-nn,:);
        Ga=Gax;G=Gx;GR=GRx;B=Bx;
    end                            %status before the last change it occured outside of the time interval

    clear Yact Yin Y0 T0 nn mk phase tc;

    %####################################################
    %###### 0 step - waiting for "equilibrium"    #######
    %####################################################

    realtime=0;
    
    while (realtime<t00)
        [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
        Yact=Y0(:,8);               %amount of NF-kBn  
        Yin=Y0(:,12);               %amount of IkBan  
        TR=Y0(:,16);                %TNF level  
        Gax=Ga;Gx=G;GRx=GR;Bx=B;
        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, calls statuschange
        
        tc=T0(mk);                  %time when the status changes 
        yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration  
        realtime=realtime+tc;
    end

    if (realtime>t00)
        nn=(realtime-t00)/dt;
        yy0=Y0(mk-nn,:);
        Ga=Gax;G=Gx;GR=GRx;B=Bx;
    end                            %status before the last change it occured ouside of the time interval

    clear Yact Yin Y0 T0 nn mk tc; 

    %#######################################
    %###### 1 step- still no TNF #######
    %#######################################
    
    realtime=0;
    ga=[Ga];g=[G];gR=[GR];                  %saves activity of IkBa A20 reporter genes 
    bb=[B];                                 %saves number of active receptors
    Y=yy0;                                  %variables where single cell run is stored
    T=zeros(1,1);
    
    while (realtime<t0)
        [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
        Yact=Y0(:,8);               %amount of NF-kBn  
        Yin=Y0(:,12);               %amount of IkBan  
        TR=Y0(:,16);                %TNF level  
        Gax=Ga;Gx=G;GRx=GR;Bx=B;
        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
        tc=T0(mk);                  %time when the status changes 
        yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
        Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
        T=[T;T0(2:mk)+realtime];     
        ga=[ga;Gax*ones(mk-1,1)];
        g=[g;Gx*ones(mk-1,1)];
        gR=[gR;GRx*ones(mk-1,1)];
        bb=[bb;Bx*ones(mk-1,1)];
        realtime=realtime+tc;
    end

    
        nn=(realtime-t0)/dt;
        x=size(Y);
        Y=Y(1:(x(1)-nn),:);
        T=T(1:(x(1)-nn));
        
        Y0(mk-nn,16)=0;        %setting TNF ON for the next step
        
        yy0=Y0(mk-nn,:);
        x1=length(ga);
        ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
        Ga=Gax;G=Gx;GR=GRx;B=Bx;
 
        
        tnfcaselist = [];
        for counter = 1:size(TNF_List,2)
            yy0(16)=10^TNF_List(counter);        %setting TNF ON for the next step
            tnfcaselist(counter,:,:) = TNFloop(t0,tw1,te1,tspan,yy0,Ga,G,GR,B,Y0,GT,AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M,T0,mk,Y,T,ga,g,gR,bb,dt);
        end        
        

%     %####################################
%     %######  2 step TNF on 1 time        #######
%     %####################################
%     
%   
%     realtime=t0;
% 
%     while (realtime<t0+tw1)
%        [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GTx=GT;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end
%     
%     
%         nn=(realtime-t0-tw1)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         Y0(mk-nn,16)=0;        %setting TNF OFF for the next step
%      
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GT=GTx;GR=GRx;B=Bx;
%  
%      
%     %########################################
%     %###### 3 step TNF washed out 1 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1;
%     
%      while (realtime<t0+tw1+te1)
%         [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%         [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%      end
%     
%    
%         nn=(realtime-t0-tw1-te1)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         %Y0(mk-nn,16)=TNF;        %setting TNF ON for the next step
%         Y0(mk-nn,16)=0;        %HT20170508-setting TNF ON for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
 
%%%%%%%%%%%%%%%%5 here
% 
%     %########################################
%     %###### 4 step TNF on for 2 time ###############
%     %########################################
%     
%     
%     realtime=t0+tw1+te1;
%     
%     while (realtime<t0+tw1+te1+tw2)
%         [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         Y0(mk-nn,16)=0;        %setting TNF OFF for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%   
%     %########################################
%     %###### 5 step TNF washed out 2 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2;
%     
%       while (realtime<t0+tw1+te1+tw2+te2)
%         [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%         [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         %Y0(mk-nn,16)=TNF;        %setting TNF ON for the next step
%         Y0(mk-nn,16)=TNF*fold*fold;        %HT20170508-setting TNF ON for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
% 
%   
%     %########################################
%     %###### 6 step TNF on for the 3 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2;
%     
%      while (realtime<t0+tw1+te1+tw2+te2+tw3)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         Y0(mk-nn,16)=0;        %setting TNF OFF for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%   
%     %########################################
%     %###### 7 step TNF washed out 3 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3;
%     
%       while (realtime<t0+tw1+te1+tw2+te2+tw3+te3)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         %Y0(mk-nn,16)=TNF;        %HT:setting TNF ON for the next step
%         Y0(mk-nn,16)=TNF*fold*fold*fold;        %setting TNF ON for the next step
% 
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%  
%      %########################################
%     %###### 8 step TNF on for the 4 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3;
%     
%      while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         Y0(mk-nn,16)=0;        %setting TNF OFF for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%   
%     %########################################
%     %###### 9 step TNF washed out 4 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4;
%     
%       while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         %Y0(mk-nn,16)=TNF;        %setting TNF ON for the next step
%         Y0(mk-nn,16)=TNF*fold*fold*fold*fold;        %setting TNF ON for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%         
%      %########################################
%     %###### 10 step TNF on for the 5 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4;
%     
%      while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         Y0(mk-nn,16)=0;        %setting TNF OFF for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%   
%     %########################################
%     %###### 11 step TNF washed out 5 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5;
%     
%       while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         %Y0(mk-nn,16)=TNF;        %setting TNF ON for the next step
%         Y0(mk-nn,16)=TNF*fold*fold*fold*fold*fold;        %setting TNF ON for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%         
%              %########################################
%     %###### 12 step TNF on for the 6 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5;
%     
%      while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5-tw6)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         Y0(mk-nn,16)=0;        %setting TNF OFF for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%   
%     %########################################
%     %###### 13 step TNF washed out 6 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6;
%     
%       while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5-tw6-te6)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         %Y0(mk-nn,16)=TNF;        %setting TNF ON for the next step
%         Y0(mk-nn,16)=TNF*fold*fold*fold*fold*fold*fold;        %setting TNF ON for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%         
%              %########################################
%     %###### 14 step TNF on for the 7 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6;
%     
%      while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5-tw6-te6-tw7)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         Y0(mk-nn,16)=0;        %setting TNF OFF for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%   
%     %########################################
%     %###### 15 step TNF washed out 7 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7;
%     
%       while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5-tw6-te6-tw7-te7)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         %Y0(mk-nn,16)=TNF;        %setting TNF ON for the next step
%         Y0(mk-nn,16)=TNF*fold*fold*fold*fold*fold*fold*fold;        %setting TNF ON for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%         
%              %########################################
%     %###### 16 step TNF on for the 8 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7;
%     
%      while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5-tw6-te6-tw7-te7-tw8)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         Y0(mk-nn,16)=0;        %setting TNF OFF for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%   
%     %########################################
%     %###### 17 step TNF washed out 8 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8;
%     
%       while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8+te8)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5-tw6-te6-tw7-te7-tw8-te8)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         %Y0(mk-nn,16)=TNF;        %setting TNF ON for the next step
%         Y0(mk-nn,16)=TNF*fold*fold*fold*fold*fold*fold*fold*fold;        %setting TNF ON for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%         
%              %########################################
%     %###### 18 step TNF on for the 9 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8+te8;
%     
%      while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8+te8+tw9)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5-tw6-te6-tw7-te7-tw8-te8-tw9)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         Y0(mk-nn,16)=0;        %setting TNF OFF for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%   
%     %########################################
%     %###### 19 step TNF washed out 9 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8+te8+tw9;
%     
%       while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8+te8+tw9+te9)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5-tw6-te6-tw7-te7-tw8-te8-tw9-te9)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         %Y0(mk-nn,16)=TNF;        %setting TNF ON for the next step
%         Y0(mk-nn,16)=TNF*fold*fold*fold*fold*fold*fold*fold*fold*fold;        %setting TNF ON for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%         
%         
%                      %########################################
%     %###### 20 step TNF on for the 10 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8+te8+tw9+te9;
%     
%      while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8+te8+tw9+te9+tw10)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5-tw6-te6-tw7-te7-tw8-te8-tw9-te9-tw10)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         Y0(mk-nn,16)=0;        %setting TNF OFF for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%   
%     %########################################
%     %###### 21 step TNF washed out 9 time ###########
%     %########################################
%     
%     
%     realtime=t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8+te8+tw9+te9+tw10;
%     
%       while (realtime<t0+tw1+te1+tw2+te2+tw3+te3+tw4+te4+tw5+te5+tw6+te6+tw7+te7+tw8+te8+tw9+te9+te9+tw10+te10)
%          [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
%         Yact=Y0(:,8);               %amount of NF-kBn  
%         Yin=Y0(:,12);               %amount of IkBan  
%         TR=Y0(:,16);                %TNF level  
%         Gax=Ga;Gx=G;GRx=GR;Bx=B;
%        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, call statuschange
%         tc=T0(mk);                  %time when the status changes 
%         yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration
%         Y=[Y;Y0(2:mk,:)];           %rows from 2 do mk, all columns
%         T=[T;T0(2:mk)+realtime];          
%         ga=[ga;Gax*ones(mk-1,1)];
%         g=[g;Gx*ones(mk-1,1)];
%         gR=[gR;GRx*ones(mk-1,1)];
%         bb=[bb;Bx*ones(mk-1,1)];
%         realtime=realtime+tc;
%     end;
%     
%    
%         nn=(realtime-t0-tw1-te1-tw2-te2-tw3-te3-tw4-te4-tw5-te5-tw6-te6-tw7-te7-tw8-te8-tw9-te9-tw10-te10)/dt;
%         x=size(Y);
%         Y=Y(1:(x(1)-nn),:);
%         T=T(1:(x(1)-nn));
%         
%         %Y0(mk-nn,16)=TNF;        %setting TNF ON for the next step
%         Y0(mk-nn,16)=0;        %setting TNF ON for the next step
%         
%         yy0=Y0(mk-nn,:);
%         x1=length(ga);
%         ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
%         Ga=Gax;G=Gx;GR=GRx;B=Bx;
%         
%         
%         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%        
       

%     YYY=YYY+Y;
%     GGa=GGa+ga;
%     GG=GG+g;
%     GGR=GGR+gR;
%     Bb=Bb+bb;
    
MM(i)=M;
NFF(i)=NF;
    
%#############################
%###### DATA FOR PLOTS #######
%#############################

    
    XXX(i,:,:,:)=Y; 
%     XB(i,:)=bb(:);   
%     XG(i,:)=g(:); 
%     XGa(i,:)=ga(:);    
%     XGR(i,:)=gR(:);  % data for all cells 
    
end                % end of the Main LOOP

% clear tspan tspan1 g ga T0 Y0 Yact Yin Y y0 yy0 tindex;
% 
% 
% 
% YYY=YYY/N;                       %average over population
% GGa=GGa/N;
% GG=GG/N;
% GGR=GGR/N;
% Bb=Bb/N;
% T=(T-t0)/60;

 
% % Data Processing : Counts Cells Responding to First and Second Pulse
% % designated for two pulses only
% 
% RecetorsAverage=sum(MM)/N
% NFaverage=sum(NFF)/N
% 
% NFKB=XXX(:,:,8)+XXX(:,:,15);
% both=0;
% onlyfirst=0;
% onlysecond=0;
% 
% both2=0;
% onlyfirst2=0;
% onlysecond2=0;
% 
% for i=1:N
%     s=0;
%     f=0;
%     s2=0;
%     f2=0;
%     aa=NFKB(i,:);
%     c1=size(aa);
%     c=round(c1(2)/2);
%     a1=aa(i:c);
%     a2=aa(c:2*c-1);
%     
% if max(a1)>meanNFkB/10  
%     f=1;
% end    
% if max(a1)>meanNFkB/5  
%     f2=1;
% end    
% if max(a2)>meanNFkB/10
%     s=1;
% end  
% if max(a2)>meanNFkB/5
%     s2=1;
% end 
%     both=both+f*s;
%     onlyfirst=onlyfirst+f-f*s;
%     onlysecond=onlysecond+s-f*s;
%     
%     both2=both2+f2*s2;
%     onlyfirst2=onlyfirst2+f2-f2*s2;
%     onlysecond2=onlysecond2+s2-f2*s2;
% end
% 
% Threshold=0.1
% 
% both
% onlyfirst
% onlysecond
% any=both+onlyfirst+onlysecond
% 
% Threshold=0.2  % Used for this study
% 
% both2
% onlyfirst2
% onlysecond2
% any2=both2+onlyfirst2+onlysecond2
% 
% %% End of data processing %%%


simulation_time=etime(clock,starttime) %simulation time seconds
save last

% %HT 2017_0508
% mcol=sqrt(N);
% mrow=sqrt(N);
% ind=1:1:N;
% X_limm=RS*60*11;
% yrightlim=TNF*fold^9;
% 
% figure(1);hold on
% for mm=1:length(ind)
% %     subplot(mrow, mcol, mm);
%     plot(T, XXX(mm,:,8))
% %     set(gca, 'xlim', [0,1000], 'ylim', [0,10*10^4], 'yscale', 'linear');
%     title('NFkB response', 'fontsize', 8);
% end
%  figure(2);clf;
% for mm=1:length(ind)
% %     subplot(mrow, mcol, mm);
%     plot(T, XXX(mm,:,16))
%     set(gca, 'xlim', [0,1000], 'ylim', [0,20], 'yscale', 'linear');
%     title('TNF conc.', 'fontsize', 8);
% end
%  
% figure(3);clf;
% for mm=1:length(ind)
% %     subplot(mrow, mcol, mm);
%     yyaxis left
%     plot(T, XXX(mm,:,8),'b')
%     set(gca, 'xlim', [0,X_limm], 'yscale', 'linear');
%     %set(gca, 'xlim', [0,X_limm], 'ylim', [0,10*10^4], 'yscale', 'linear');
%     ylabel('Nuc. P65')
%     hold on
%     yyaxis right
%     plot(T, XXX(mm,:,16),'r')
%     set(gca, 'xlim', [0,X_limm], 'ylim', [0,yrightlim], 'yscale', 'linear');
%     ylabel('Ligand Conc.')
% end
%AllCellPlotting
>>>>>>> e409751e531eb57f3e9bdec15b37d1cdd56cfeee
