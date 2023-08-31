close all;
clc;
clear all;

%% Clarification: this is another test exmaple we have done for your reference. it's sad we lost data due to weather. We will try to update but cannot guarantee. 

global  A B1 B2 B3 u1 u2 u3 u1e u2e u3e Kl1 Kl2 Kl3 R1 R2 R3 S1 S2 S3

A=[-3 -1;
    1 -3];
% B1=[-1;1];
% B2=[1;-1];
% B3=[1;1];

B1=[0  1]';
B2=[1  0]';
B3=[1  1]';


x=[0.01;-0.01];
xe=[0.01;-0.01];



i=1;
ii=0;



K1=[0.2 0.9];
K2=[0.3 0.05];
K3=[0.3 0.1];

K1i=K1;
K2i=K2;
K3i=K3;


K1e=[0.5 1];
K2e=[0.5 0.125];
K3e=[0.4 0.4];

R1=1;
R2=2;
R3=3;

S1i=eye(2);
S2i=eye(2);
S3i=eye(2);


P1i=zeros(2,2);
P2i=zeros(2,2);
P3i=zeros(2,2);

dKl1=[];
dKl2=[];
dKl3=[];


dKl1i=1;


T=0.001;

for j=1:1

for k=1:2000
    
    Kl1last=K1i(:,:,end);
    Kl2last=K2i(:,:,end);
    Kl3last=K3i(:,:,end);
    
    
    e1=0.0001*rand(1);
    e2=0.0001*rand(1);
    e3=0.0001*rand(1);
    
    
    e4=0.0000*rand(1);
    e5=0.0000*rand(1);
    e6=0.0000*rand(1);
    
    u1=-K1*x(:,k)+e1;
    u2=-K2*x(:,k)+e2;
    u3=-K3*x(:,k)+e3;
    
    Kl1=K1i(:,:,i);
    Kl2=K2i(:,:,i);
    Kl3=K3i(:,:,i);
    
    S1=S1i(:,:,i);
    S2=S2i(:,:,i);
    S3=S3i(:,:,i);
    
    
    tspanxd=[0 T];
    X(:,k)=[x(:,k)',zeros(1,10)]';
    [tl,dl]= ode45(@learner,tspanxd,X(:,k));
    x(:,k+1)=[dl(length(tl),1); dl(length(tl),2)];
    
    
    u1e=-K1e*xe(:,k)+e4;
    u2e=-K2e*xe(:,k)+e5;
    u3e=-K3e*xe(:,k)+e6;
    
    Xe(:,k)=[xe(:,k)',zeros(1,13)]';
    [te,de]= ode45(@expert,tspanxd,Xe(:,k));
    xe(:,k+1)=[de(length(te),1); de(length(te),2)];   
    
    ii=ii+1;
    dxx(ii,:)=[x(1,k+1)^2,2*x(1,k+1)*x(2,k+1),x(2,k+1)^2]-[x(1,k)^2,2*x(1,k)*x(2,k),x(2,k)^2]+[xe(1,k+1)^2,2*xe(1,k+1)*xe(2,k+1),xe(2,k+1)^2]-[xe(1,k)^2,2*xe(1,k)*xe(2,k),xe(2,k)^2];
    
    lxu1(ii,:)=[-dl(length(tl),3)-de(length(te),3) -dl(length(tl),4)-de(length(te),4)];
    lxu2(ii,:)=[-dl(length(tl),5)-de(length(te),5) -dl(length(tl),6)-de(length(te),6)];
    lxu3(ii,:)=[-dl(length(tl),7)-de(length(te),7) -dl(length(tl),8)-de(length(te),8)];
    
    lxx1(ii)=-dl(length(tl),9); 
    lxx2(ii)=-dl(length(tl),10);
    lxx3(ii)=-dl(length(tl),11);
    
    
    r1(ii)= -de(length(te),9)-de(length(te),12)-de(length(te),13)-dl(length(tl),12); 
         
    r2(ii)=-de(length(te),10)-de(length(te),12)-de(length(te),14)-dl(length(tl),12);
         
    r3(ii)=-de(length(te),11)-de(length(te),12)-de(length(te),15)-dl(length(tl),12);
         
     
    fai(ii,:)=[dxx(ii,:),lxu1(ii,:),lxu2(ii,:),lxu3(ii,:),lxx1(ii),lxx2(ii),lxx3(ii)];     

    
    
    a=rank(fai'*fai)
    
    if(a==12 && dKl1i>0.001)
        
        
        p1=(fai'*fai)\fai'*r1';
        p2=(fai'*fai)\fai'*r2';
        p3=(fai'*fai)\fai'*r3';
        
       % player 1
        P1(:,:,i)=[p1(1) p1(2);
                   p1(2) p1(3)];
        K1i(:,:,i+1)=inv(R1)*[p1(4) p1(5)]'; 
        
        S1i(:,:,i+1)=[p1(10) p1(11);
                      p1(11) p1(12)];
        
        % player 2       
        P2(:,:,i)=[p2(1) p2(2);
                   p2(2) p2(3)];
        K2i(:,:,i+1)=inv(R2)*[p2(6) p2(7)]'; 
        
        S2i(:,:,i+1)=[p2(10) p2(11);
                      p2(11) p2(12)];       
               
        % player 3       
               
        P3(:,:,i)=[p3(1) p3(2);
                   p3(2) p3(3)];
        K3i(:,:,i+1)=inv(R3)*[p3(8) p3(9)]'; 
        
        S3i(:,:,i+1)=[p3(10) p3(11);
                      p3(11) p3(12)];    
               
               
        dKl1i=norm(K1i(:,:,i+1)-K1i(:,:,i));       
        dKl1=[dKl1,norm(K1i(:,:,i+1)-Kl1last)];
        dKl2=[dKl2,norm(K2i(:,:,i+1)-Kl2last)];
        dKl3=[dKl3,norm(K3i(:,:,i+1)-Kl3last)]; %testing convergence of K





        
            ii=0;
            dxx=[];
            lxu1=[];
            lxu2=[];
            lxu3=[];
            r1=[];
            r2=[];
            r3=[];
            fai=[];
            a=0;
        
        
        
       i=i+1;
    
    end
    
     if(dKl1i(end)<=0.001)
            P1i(:,:,j)=P1(:,:,end);
            P2i(:,:,j)=P2(:,:,end);
            P3i(:,:,j)=P3(:,:,end);
            Kl1i(:,:,j)=K1i(:,end);
            Kl2i(:,:,j)=K2i(:,end);
            Kl3i(:,:,j)=K3i(:,end);
%             S1i(:,:,j)=S1(:,:,end);
%             S2i(:,:,j)=S2(:,:,end);
%             S3i(:,:,j)=S3(:,:,end);
          

           break;
        end
    
    
    
end


end




