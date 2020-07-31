<<<<<<< HEAD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                           %
    %         Function includes system of ODEs describing       %
    %         NF-kB regulatory pathway.                         %
    %         Substrates are coded as follows:                  %
    %                                                           %  
    %         y(1)   IKKKa active                               %
    %         y(2)   IKKn   neutral                             %    
    %         y(3)   IKKa   active                              %
    %         y(4)   IKKi   inactive                            %
    %         y(5)   phospho-IkBa cytoplasmic                   %
    %         y(6)   phospho-IkBa|NFkB cytoplasmic              %  
    %         y(7)   NFkB  cytoplasmic                          %
    %         y(8)   NFkBn  nuclear                             %
    %         y(9)   A20                                        %
    %         y(10)  A20t                                       %
    %         y(11)  IkBa                                       %
    %         y(12)  IkBan                                      %
    %         y(13)  IkBat                                      %      
    %         y(14)  (IkBa|NFkB) cytoplasmic                    %  
    %         y(15)  (IkBan|NFkBn) nuclear                      %   
    %         y(16)  extracellular TNF                          %   
    %                                                           %
    %         y(17)  mRNA reporter                              %
    %                                                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function dy=Model(t,y,Ga,G,GR,B)

 [NF0,NF1,NF2,M0,M1,M2,k4,ka20,AB,kv,q1,q2,c1,c3,c4,c5,k1,k2,k3,a1,a2,a3,c1a,c5a,c6a,i1,i1a,e1a,e2a,dt,tp,KN,KNN,ka,ki,kb,kf,Tdeg,q1r,q2r,q2rr,c1r,c1rr,c3r]=Parameters;
  
 %###############################################################
 
 dy=zeros(19,1);
 
 dy(1)=ka*B*(KN-y(1))* ka20/(ka20+y(9))-ki*y(1);                            %active IKKK kinase 
 dy(2)=-y(1)^2*k1*y(2)+k4*(KNN-y(2)-y(3)-y(4));                                 %neutral IKK   
 dy(3)=y(1)^2*k1*y(2)-k3*y(3)*(k2+y(9))/k2;                                     %free active IKK                                                                                    
 dy(4)=k3*y(3)*(k2+y(9))/k2-k4*y(4);                                            %inactive IKK   
 dy(5)=a2*y(3)*y(11)-tp*y(5);                                                   %Phospo-IkBa cytoplasmic 
 dy(6)=a3*y(3)*y(14)-tp*y(6);                                                   %cytoplasmic (phospho-IkBa|NF-kB) 
 dy(7)=c6a*y(14)-a1*y(7)*y(11)+tp*y(6)-i1*y(7);                                 %free cytoplasmic NFkB
 dy(8)=i1*y(7)-a1*kv*y(12)*y(8);                                                %free nuclear NFkB
 dy(9)=c4*y(10)-c5*y(9);                                                        %cytoplasmic A20
 dy(10)=c1*G-c3*y(10);                                                          %A20 transcript
 dy(11)=-a2*y(3)*y(11)-a1*y(11)*y(7)+c4*y(13)-c5a*y(11)-i1a*y(11)+e1a*y(12);    %free cytoplasmic IkBa
 dy(12)=-a1*kv*y(12)*y(8)+i1a*y(11)-e1a*y(12);                                  %free nuclear IkBan
 dy(13)=c1a*Ga-c3*y(13);                                                        %IkBa transcript
 dy(14)=a1*y(11)*y(7)-c6a*y(14)-a3*y(3)*y(14)+e2a*y(15);                        %cytoplasmic (IkBa|NFkB) complex
 dy(15)=a1*kv*y(12)*y(8)-e2a*y(15);                                             %nuclear (IkBa|NFkB) complex
 dy(16)=Tdeg*y(16);                                                            %extracellular TNF
 dy(17)=c1rr+c1r*GR-c3r*y(17);                                                  %Reporter transcript
 
 
=======
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                           %
    %         Function includes system of ODEs describing       %
    %         NF-kB regulatory pathway.                         %
    %         Substrates are coded as follows:                  %
    %                                                           %  
    %         y(1)   IKKKa active                               %
    %         y(2)   IKKn   neutral                             %    
    %         y(3)   IKKa   active                              %
    %         y(4)   IKKi   inactive                            %
    %         y(5)   phospho-IkBa cytoplasmic                   %
    %         y(6)   phospho-IkBa|NFkB cytoplasmic              %  
    %         y(7)   NFkB  cytoplasmic                          %
    %         y(8)   NFkBn  nuclear                             %
    %         y(9)   A20                                        %
    %         y(10)  A20t                                       %
    %         y(11)  IkBa                                       %
    %         y(12)  IkBan                                      %
    %         y(13)  IkBat                                      %      
    %         y(14)  (IkBa|NFkB) cytoplasmic                    %  
    %         y(15)  (IkBan|NFkBn) nuclear                      %   
    %         y(16)  extracellular TNF                          %   
    %                                                           %
    %         y(17)  mRNA reporter                              %
    %                                                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function dy=Model(t,y,Ga,G,GR,B)

 [NF0,NF1,NF2,M0,M1,M2,k4,ka20,AB,kv,q1,q2,c1,c3,c4,c5,k1,k2,k3,a1,a2,a3,c1a,c5a,c6a,i1,i1a,e1a,e2a,dt,tp,KN,KNN,ka,ki,kb,kf,Tdeg,q1r,q2r,q2rr,c1r,c1rr,c3r]=Parameters;
  
 %###############################################################
 
 dy=zeros(19,1);
 
 dy(1)=ka*B*(KN-y(1))* ka20/(ka20+y(9))-ki*y(1);                            %active IKKK kinase 
 dy(2)=-y(1)^2*k1*y(2)+k4*(KNN-y(2)-y(3)-y(4));                                 %neutral IKK   
 dy(3)=y(1)^2*k1*y(2)-k3*y(3)*(k2+y(9))/k2;                                     %free active IKK                                                                                    
 dy(4)=k3*y(3)*(k2+y(9))/k2-k4*y(4);                                            %inactive IKK   
 dy(5)=a2*y(3)*y(11)-tp*y(5);                                                   %Phospo-IkBa cytoplasmic 
 dy(6)=a3*y(3)*y(14)-tp*y(6);                                                   %cytoplasmic (phospho-IkBa|NF-kB) 
 dy(7)=c6a*y(14)-a1*y(7)*y(11)+tp*y(6)-i1*y(7);                                 %free cytoplasmic NFkB
 dy(8)=i1*y(7)-a1*kv*y(12)*y(8);                                                %free nuclear NFkB
 dy(9)=c4*y(10)-c5*y(9);                                                        %cytoplasmic A20
 dy(10)=c1*G-c3*y(10);                                                          %A20 transcript
 dy(11)=-a2*y(3)*y(11)-a1*y(11)*y(7)+c4*y(13)-c5a*y(11)-i1a*y(11)+e1a*y(12);    %free cytoplasmic IkBa
 dy(12)=-a1*kv*y(12)*y(8)+i1a*y(11)-e1a*y(12);                                  %free nuclear IkBan
 dy(13)=c1a*Ga-c3*y(13);                                                        %IkBa transcript
 dy(14)=a1*y(11)*y(7)-c6a*y(14)-a3*y(3)*y(14)+e2a*y(15);                        %cytoplasmic (IkBa|NFkB) complex
 dy(15)=a1*kv*y(12)*y(8)-e2a*y(15);                                             %nuclear (IkBa|NFkB) complex
 dy(16)=Tdeg*y(16);                                                            %extracellular TNF
 dy(17)=c1rr+c1r*GR-c3r*y(17);                                                  %Reporter transcript
 
 
>>>>>>> e409751e531eb57f3e9bdec15b37d1cdd56cfeee
  