close all;
clc;
clear all;



% global A B1 B2  S1 S2   Q1 Q2 R1 R2;

A=[-3 -1;
   1 -3];
B1=[0;1];
B2=[1;0];
B3=[1;0];



% A=[-1.01887 0.90506 -0.00215;
%     0.82225 -1.07741 -0.1755;
%     0        0          -1   ];
% B1=[0 0  2]';
% B2=[0 0  2]';
% B3=[0 0  3]';
%% expert parameters 
% Pe1= [1  0 0 ;
%       0  1 0 ;
%       0 0  1];%[0.5 0; 0  0.5];
% Pe2= [1  0 0;
%       0  2 0;
%       0  0  1];%[0.5 0; 0  1]
%   Pe3= [2  0 0;
%         0  2 0;
%         0  0  2];


Pe1= [1  0;
      0  1];
Pe2= [1  0;
      0  2];
Pe3= [0.5  0;
      0  1];


R1e=1;
R2e=1;
R3e=1;


S1e=-((A-B1*inv(R1e)*B1'*Pe1-B2*inv(R2e)*B2'*Pe2-B3*inv(R3e)*B3'*Pe3)'*Pe1+Pe1*(A-B1*inv(R1e)*B1'*Pe1-B2*inv(R2e)*B2'*Pe2-B3*inv(R3e)*B3'*Pe3)+Pe1*B1*inv(R1e)*B1'*Pe1+Pe2*B2*inv(R2e)*B2'*Pe2+Pe3*B3*inv(R3e)*B3'*Pe3);
S2e=-((A-B1*inv(R1e)*B1'*Pe1-B2*inv(R2e)*B2'*Pe2-B3*inv(R3e)*B3'*Pe3)'*Pe2+Pe2*(A-B1*inv(R1e)*B1'*Pe1-B2*inv(R2e)*B2'*Pe2-B3*inv(R3e)*B3'*Pe3)+Pe1*B1*inv(R1e)*B1'*Pe1+Pe2*B2*inv(R2e)*B2'*Pe2+Pe3*B3*inv(R3e)*B3'*Pe3);
S3e=-((A-B1*inv(R1e)*B1'*Pe1-B2*inv(R2e)*B2'*Pe2-B3*inv(R3e)*B3'*Pe3)'*Pe3+Pe3*(A-B1*inv(R1e)*B1'*Pe1-B2*inv(R2e)*B2'*Pe2-B3*inv(R3e)*B3'*Pe3)+Pe1*B1*inv(R1e)*B1'*Pe1+Pe2*B2*inv(R2e)*B2'*Pe2+Pe3*B3*inv(R3e)*B3'*Pe3);



Q1e=[A'*Pe1+Pe1*A+Pe1+S1e Pe1*B1 Pe1*B2  Pe1*B3;
       B1'*Pe1             R1e      0      0;
       B2'*Pe1              0       R2e    0;
       B3'*Pe1              0        0     R3e];
   
 
Q2e=[A'*Pe2+Pe2*A+Pe2+S2e Pe2*B1 Pe2*B2  Pe2*B3;
       B1'*Pe2             R1e      0      0;
       B2'*Pe2              0       R2e    0
       B3'*Pe2              0       0      R3e];

Q3e=[A'*Pe3+Pe3*A+Pe3+S3e Pe3*B1 Pe3*B2  Pe3*B3;
       B1'*Pe3             R1e      0      0;
       B2'*Pe3              0       R2e    0
       B3'*Pe3              0       0      R3e];

Q1e=[Q1e(1,1) Q1e(1,2) Q1e(1,3) Q1e(1,4) Q1e(1,5) Q1e(2,2) Q1e(2,3) Q1e(2,4) Q1e(2,5) Q1e(3,3) Q1e(3,4) Q1e(3,5) Q1e(4,4) Q1e(4,5) Q1e(5,5)];
Q2e=[Q2e(1,1) Q2e(1,2) Q2e(1,3) Q2e(1,4) Q2e(1,5) Q2e(2,2) Q2e(2,3) Q2e(2,4) Q2e(2,5) Q2e(3,3) Q2e(3,4) Q2e(3,5) Q2e(4,4) Q2e(4,5) Q2e(5,5)];
Q3e=[Q3e(1,1) Q3e(1,2) Q3e(1,3) Q3e(1,4) Q3e(1,5) Q3e(2,2) Q3e(2,3) Q3e(2,4) Q3e(2,5) Q3e(3,3) Q3e(3,4) Q3e(3,5) Q3e(4,4) Q3e(4,5) Q3e(5,5)];



K1e=inv(Q1e(10))*[Q1e(3), Q1e(7)];
K2e=inv(Q2e(13))*[Q2e(4), Q2e(8)];
K3e=inv(Q3e(15))*[Q3e(5), Q3e(9)];

% K1e=inv(Q1e(8))*[Q1e(3), Q1e(6)];
% K2e=inv(Q2e(10))*[Q2e(4), Q2e(7)];
% K2e=[0 2];



%  S1=[3,     0;
%      0,   3.25];
%  S2=[6,    0;
%      0,  10.75];
%   S1=[1,0;
%      0, 1];
%  S2=[1,0;
%     0, 1];  



  S1=eye(2);
 S2=eye(2);
  S3=eye(2);
   



R1=1*eye(1);
R2=1*eye(1);
R3=1*eye(1);


K1=[0.01 0.5];
K2=[0.01 1.1];
K3=[0.01 2.1];


% Q1=[3.25    0;
%     0  3.25];
%  
%  Q2=[3.25     0;
%      0   7.25];


 iter=2000;



for i=1:iter
    
   
    
    %% use expert information to update kernel matrix Q1, Q2
     
       Abar(:,:,i)=A-B1*K1(:,:,i)-B2*K2(:,:,i)-B3*K3(:,:,i);
       Q1(:,:,i)=K1(:,:,i)'*R1*K1(:,:,i)+K2(:,:,i)'*R2*K2(:,:,i)+K3(:,:,i)'*R3*K3(:,:,i)+(K1(:,:,i)-K1e)'*R1*(K1(:,:,i)-K1e)+S1(:,:,i);
       Q2(:,:,i)=K1(:,:,i)'*R1*K1(:,:,i)+K2(:,:,i)'*R2*K2(:,:,i)+K3(:,:,i)'*R3*K3(:,:,i)+(K2(:,:,i)-K2e)'*R2*(K2(:,:,i)-K2e)+S2(:,:,i);       
       Q3(:,:,i)=K1(:,:,i)'*R1*K1(:,:,i)+K2(:,:,i)'*R2*K2(:,:,i)+K3(:,:,i)'*R3*K3(:,:,i)+(K3(:,:,i)-K3e)'*R3*(K3(:,:,i)-K3e)+S3(:,:,i);       
       
       P1(:,:,i)=lyap(Abar(:,:,i)',Q1(:,:,i));
       P2(:,:,i)=lyap(Abar(:,:,i)',Q2(:,:,i));
       P3(:,:,i)=lyap(Abar(:,:,i)',Q3(:,:,i));
              
       S1(:,:,i+1)=-(A-B1*K1(:,:,i)-B2*K2(:,:,i)-B3*K3(:,:,i))'*P1(:,:,i)-P1(:,:,i)*(A-B1*K1(:,:,i)-B2*K2(:,:,i)-B3*K3(:,:,i))-K1(:,:,i)'*R1*K1(:,:,i)-K2(:,:,i)'*R2*K2(:,:,i)-K3(:,:,i)'*R3*K3(:,:,i);
       S2(:,:,i+1)=-(A-B1*K1(:,:,i)-B2*K2(:,:,i)-B3*K3(:,:,i))'*P2(:,:,i)-P2(:,:,i)*(A-B1*K1(:,:,i)-B2*K2(:,:,i)-B3*K3(:,:,i))-K1(:,:,i)'*R1*K1(:,:,i)-K2(:,:,i)'*R2*K2(:,:,i)-K3(:,:,i)'*R3*K3(:,:,i);
       S3(:,:,i+1)=-(A-B1*K1(:,:,i)-B2*K2(:,:,i)-B3*K3(:,:,i))'*P3(:,:,i)-P3(:,:,i)*(A-B1*K1(:,:,i)-B2*K2(:,:,i)-B3*K3(:,:,i))-K1(:,:,i)'*R1*K1(:,:,i)-K2(:,:,i)'*R2*K2(:,:,i)-K3(:,:,i)'*R3*K3(:,:,i);
            
       K1(:,:,i+1)=inv(R1)*B1'*P1(:,:,i);
       K2(:,:,i+1)=inv(R2)*B2'*P2(:,:,i);
       K3(:,:,i+1)=inv(R3)*B3'*P3(:,:,i);
       
       
       p1(:,i)=norm(P1(:,:,i));
       p2(:,i)=norm(P2(:,:,i));
       p3(:,i)=norm(P3(:,:,i));
       
       s1(:,i)=norm(S1(:,:,i));
       s2(:,i)=norm(S2(:,:,i));
       s3(:,i)=norm(S3(:,:,i));
       
       k1(:,i)=norm(K1(:,:,i)-K1e);
       k2(:,i)=norm(K2(:,:,i)-K2e);
       k3(:,i)=norm(K3(:,:,i)-K3e);
       
       
       
       
              
end 
       
     figure  
plot(p1,'LineWidth',2);
xlabel('Time iteration $h$');ylabel(' $\Vert P_1^h \Vert$','LineWidth',2);
set(get(gca,'XLabel'),'FontSize',15);set(get(gca,'YLabel'),'FontSize',15);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');

   
     figure  
plot(p2,'LineWidth',2);
xlabel('Time iteration $h$');ylabel(' $\Vert P_2^h \Vert$','LineWidth',2);
set(get(gca,'XLabel'),'FontSize',15);set(get(gca,'YLabel'),'FontSize',15);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');

     figure  
plot(p3,'LineWidth',2);
xlabel('Time iteration $h$');ylabel(' $\Vert P_3^h \Vert$','LineWidth',2);
set(get(gca,'XLabel'),'FontSize',15);set(get(gca,'YLabel'),'FontSize',15);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');

   
     figure  
plot(k1,'LineWidth',2); 
xlabel('Time iteration $h$');ylabel(' $\Vert K_1^h-K_{1e} \Vert$','LineWidth',2);
set(get(gca,'XLabel'),'FontSize',15);set(get(gca,'YLabel'),'FontSize',15);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');

   
     figure  
plot(k2,'LineWidth',2);
xlabel('Time iteration $h$');ylabel(' $\Vert K_2^h-K_{2e} \Vert$','LineWidth',2);
set(get(gca,'XLabel'),'FontSize',15);set(get(gca,'YLabel'),'FontSize',15);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');

    figure  
plot(k3,'LineWidth',2);
xlabel('Time iteration $h$');ylabel(' $\Vert K_3^h-K_{3e} \Vert$','LineWidth',2);
set(get(gca,'XLabel'),'FontSize',15);set(get(gca,'YLabel'),'FontSize',15);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');

   
     figure  
plot(s1,'LineWidth',2); 
xlabel('Time iteration $h$');ylabel(' $\Vert S_1^h \Vert$','LineWidth',2);
set(get(gca,'XLabel'),'FontSize',15);set(get(gca,'YLabel'),'FontSize',15);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');


   
     figure  
plot(s2,'LineWidth',2); 
xlabel('Time iteration $h$');ylabel(' $\Vert S_2^h \Vert$','LineWidth',2);
set(get(gca,'XLabel'),'FontSize',15);set(get(gca,'YLabel'),'FontSize',15);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');
  
     figure  
plot(s3,'LineWidth',2); 
xlabel('Time iteration $h$');ylabel(' $\Vert S_3^h \Vert$','LineWidth',2);
set(get(gca,'XLabel'),'FontSize',15);set(get(gca,'YLabel'),'FontSize',15);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');
