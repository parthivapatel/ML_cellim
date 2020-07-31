<<<<<<< HEAD
function [Y]=TNFloop(t0,tw1,te1,tspan,yy0,Ga,G,GR,B,Y0,GT,AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M,T0,mk,Y,T,ga,g,gR,bb,dt)  


realtime=t0;

while (realtime<t0+tw1)
   [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
    Yact=Y0(:,8);               %amount of NF-kBn  
    Yin=Y0(:,12);               %amount of IkBan  
    TR=Y0(:,16);                %TNF level  
    Gax=Ga;Gx=G;GTx=GT;GRx=GR;Bx=B;
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


nn=(realtime-t0-tw1)/dt;
x=size(Y);
Y=Y(1:(x(1)-nn),:);
T=T(1:(x(1)-nn));

Y0(mk-nn,16)=0;        %setting TNF OFF for the next step

yy0=Y0(mk-nn,:);

x1=length(ga);
ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
Ga=Gax;G=Gx;GT=GTx;GR=GRx;B=Bx;


realtime=t0+tw1;
    
 while (realtime<t0+tw1+te1)
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


nn=(realtime-t0-tw1-te1)/dt;
x=size(Y);
Y=Y(1:(x(1)-nn),:);
T=T(1:(x(1)-nn));

%Y0(mk-nn,16)=TNF;        %setting TNF ON for the next step
Y0(mk-nn,16)=0;        %HT20170508-setting TNF ON for the next step

yy0=Y0(mk-nn,:);
x1=length(ga);
ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
Ga=Gax;G=Gx;GR=GRx;B=Bx;

=======
function [Y]=TNFloop(t0,tw1,te1,tspan,yy0,Ga,G,GR,B,Y0,GT,AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M,T0,mk,Y,T,ga,g,gR,bb,dt)  


realtime=t0;

while (realtime<t0+tw1)
   [T0,Y0]=ode23tb(@Model,tspan,yy0,[],Ga,G,GR,B); 
    Yact=Y0(:,8);               %amount of NF-kBn  
    Yin=Y0(:,12);               %amount of IkBan  
    TR=Y0(:,16);                %TNF level  
    Gax=Ga;Gx=G;GTx=GT;GRx=GR;Bx=B;
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


nn=(realtime-t0-tw1)/dt;
x=size(Y);
Y=Y(1:(x(1)-nn),:);
T=T(1:(x(1)-nn));

Y0(mk-nn,16)=0;        %setting TNF OFF for the next step

yy0=Y0(mk-nn,:);

x1=length(ga);
ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
Ga=Gax;G=Gx;GT=GTx;GR=GRx;B=Bx;


realtime=t0+tw1;
    
 while (realtime<t0+tw1+te1)
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


nn=(realtime-t0-tw1-te1)/dt;
x=size(Y);
Y=Y(1:(x(1)-nn),:);
T=T(1:(x(1)-nn));

%Y0(mk-nn,16)=TNF;        %setting TNF ON for the next step
Y0(mk-nn,16)=0;        %HT20170508-setting TNF ON for the next step

yy0=Y0(mk-nn,:);
x1=length(ga);
ga=ga(1:x1-nn);g=g(1:x1-nn);gR=gR(1:x1-nn);bb=bb(1:x1-nn);
Ga=Gax;G=Gx;GR=GRx;B=Bx;

>>>>>>> e409751e531eb57f3e9bdec15b37d1cdd56cfeee
end