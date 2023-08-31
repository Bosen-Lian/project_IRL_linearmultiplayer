function xdot=expert(t,x)
global A;
global B1 B2 B3;
global u1e u2e u3e;
global Kl1 Kl2 Kl3;
global R1 R2 R3;
global S1 S2 S3;


x=[x(1);x(2)];


xdot=[A*x+B1*u1e+B2*u2e+B3*u3e       
      2*(u1e+Kl1*x)*x
      2*(u2e+Kl2*x)*x
      2*(u3e+Kl3*x)*x
       x'*S1*x % 9-th
       x'*S2*x
       x'*S3*x
      (Kl1*x)'*R1*Kl1*x+(Kl2*x)'*R2*Kl2*x+(Kl3*x)'*R3*Kl3*x %12-th
      (u1e+Kl1*x)'*R1*(u1e+Kl1*x)
      (u2e+Kl2*x)'*R2*(u2e+Kl2*x)
      (u3e+Kl3*x)'*R3*(u3e+Kl3*x)];
   
     

end