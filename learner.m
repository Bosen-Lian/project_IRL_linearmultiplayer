function xdot=learner(t,x)
global A;
global B1 B2 B3;
global u1 u2 u3;
global Kl1 Kl2 Kl3;
global R1 R2 R3;


x=[x(1);x(2)];



% x'*Kl1'*R1*Kl1*x

xdot=[A*x+B1*u1+B2*u2+B3*u3       
      2*(u1+Kl1*x)*x
      2*(u2+Kl2*x)*x
      2*(u3+Kl3*x)*x
      x(1)*x(1) %9-th
      2*x(1)*x(2)
      x(2)*x(2)
      (Kl1*x)'*R1*Kl1*x+(Kl2*x)'*R2*Kl2*x+(Kl3*x)'*R3*Kl3*x];
   
     

end