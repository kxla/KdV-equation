function Ex1
clc;tic;
baita=6;geima=1;
XL=-40;XR=60;h=0.05;tao=0.001;Tfinal=1;
M=ceil((XR-XL)/h);M1=M+1;M3=M-3;
T=ceil(Tfinal/tao);T1=T+1;
for j=1:M1
    x(j)=XL+(j-1)*h;
end
for n=1:T1
    t(n)=(n-1)*tao;
end
% Exact Solution
for j=1:M1
    uexact(j)=12*(3+4*cosh(2*x(j)-8*t(T1))+cosh(4*x(j)-64*t(T1)))/((3*cosh(x(j)-28*t(T1))+cosh(3*x(j)-36*t(T1)))^2);
end     
% Numerical Solution
% Initial Conditions
for j=1:M1
    u(j,1)=12*(3+4*cosh(2*x(j))+cosh(4*x(j)))/((3*cosh(x(j))+cosh(3*x(j)))^2);
end
% Start the iterations
%Boundary conditions
for n=1:T1
    u(1,n)=12*(3+4*cosh(2*x(1)-8*t(n))+cosh(4*x(1)-64*t(n)))/((3*cosh(x(1)-28*t(n))+cosh(3*x(1)-36*t(n)))^2);
    u(2,n)=12*(3+4*cosh(2*x(2)-8*t(n))+cosh(4*x(2)-64*t(n)))/((3*cosh(x(2)-28*t(n))+cosh(3*x(2)-36*t(n)))^2);
    u(M,n)=12*(3+4*cosh(2*x(M)-8*t(n))+cosh(4*x(M)-64*t(n)))/((3*cosh(x(M)-28*t(n))+cosh(3*x(M)-36*t(n)))^2);
    u(M1,n)=12*(3+4*cosh(2*x(M1)-8*t(n))+cosh(4*x(M1)-64*t(n)))/((3*cosh(x(M1)-28*t(n))+cosh(3*x(M1)-36*t(n)))^2);
end
for j=1:M3
    for n=1:T1
        v(j,n)=0;
    end
end
%step 2
%solve U2
A=zeros(M3,M3);B=zeros(M3,M3);
A0=-geima*tao;
B0=2*geima*tao+h^3;
C0=2*h^3;
D0=-2*geima*tao+h^3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:M3
    E0(j)=-A0+baita*tao*h^2/12*u(j+2-2,1);
    F0(j)=D0+5*baita*tao*h^2/6*u(j+2-1,1);
    G0(j)=B0-5*baita*tao*h^2/6*u(j+2+1,1);
    H0(j)=A0-baita*tao*h^2/12*u(j+2+2,1);
end
for i=3:M3
    A(i,i-2)=A0;B(i,i-2)=E0(i);%u(j-2)
end
for i=2:M3
    A(i,i-1)=B0;B(i,i-1)=F0(i);%u(j-1)
end
for i=1:M3
    A(i,i)=C0;B(i,i)=C0;%u(j)
end
for i=1:(M3-1)
    A(i,i+1)=D0;B(i,i+1)=G0(i);%u(j+1)
end
for i=1:(M3-2)
    A(i,i+2)=-A0;B(i,i+2)=H0(i);%u(j+2)
end
for i=1:M3
    s(i)=0;
    for j=1:M3
        s(i)=s(i)+B(i,j)*u(j+2,1);
    end
    F(i)=s(i);
end
Y0=F';
% Solve a pentadiagonal system Au=F
% Extract bands
d=diag(A);
a=[diag(A,1);0];
b=[diag(A,2);0;0];
c=[0;diag(A,-1)];
e=[0;0;diag(A,-2)];
% Solve a pentadiagonal system Au=F
% Extract bands
N=length(Y0);        
alpha=zeros(N-1,1);beta=zeros(N-2,1);z=zeros(N,1);
gama=zeros(N,1);mv=zeros(N,1);
gama(1)=0;           
% Factor A=LR
mv(1)=d(1);alpha(1)=a(1)/mv(1);
beta(1)=b(1)/mv(1);z(1)=Y0(1)/mv(1);
gama(2)=c(2);
mv(2)=d(2)-alpha(1)*gama(2);
alpha(2)=(a(2)-beta(1)*gama(2))/mv(2);
beta(2)=b(2)/mv(2);
z(2)=(Y0(2)-z(1)*gama(2))/mv(2);  
for i=3:N-2
    gama(i)=c(i)-alpha(i-2)*e(i);
    mv(i)=d(i)-beta(i-2)*e(i)-alpha(i-1)*gama(i);
    alpha(i)=(a(i)-beta(i-1)*gama(i))/mv(i);
    beta(i)=b(i)/mv(i);
    z(i)=(Y0(i)-z(i-2)*e(i)-z(i-1)*gama(i))/mv(i);
end    
gama(N-1)=c(N-1)-alpha(N-3)*e(N-1);
mv(N-1)=d(N-1)-beta(N-3)*e(N-1)-alpha(N-2)*gama(N-1);
alpha(N-1)=(a(N-1)-beta(N-2)*gama(N-1))/mv(N-1);
gama(N)=c(N)-alpha(N-2)*e(N);
mv(N)=d(N)-beta(N-2)*e(N)-alpha(N-1)*gama(N);
z(N-1)=(Y0(N-1)-z(N-3)*e(N-1)-z(N-2)*gama(N-1))/mv(N-1);
z(N)=(Y0(N)-z(N-2)*e(N)-z(N-1)*gama(N))/mv(N);  
% Back substitution
v(N,2)=z(N);v(N-1,2)=z(N-1)-alpha(N-1)*v(N,2);
for k=N-2:-1:1
    v(k,2)=z(k)-alpha(k)*v(k+1,2)-beta(k)*v(k+2,2);    
end
for j=3:(M-1)
    u(j,2)=v(j-2,2);
end
%step 3
C1=2*h^3;
for n=3:T1
    for j=1:M3
        A1(j)=-tao*(2*geima+h^2/6*0.5*baita*u(j+2-2,n-1));
        B1(j)=4*tao*geima+h^3-5*tao*baita*h^2/6*u(j+2-1,n-1);
        D1(j)=-4*tao*geima+h^3+5*tao*baita*h^2/6*u(j+2+1,n-1);
        E1(j)=tao*(2*geima+h^2/6*0.5*baita*u(j+2+2,n-1));
        F1(j)=-4*tao*geima+h^3+5*tao*baita*h^2/6*u(j+2-1,n-1);
        G1(j)=4*tao*geima+h^3-5*tao*baita*h^2/6*u(j+2+1,n-1);
    end
    for i=3:M3
        A(i,i-2)=A1(i);B(i,i-2)=-A1(i);
    end 
    for i=2:M3
        A(i,i-1)=B1(i);B(i,i-1)=F1(i);
    end    
    for i=1:M3
        A(i,i)=C1;B(i,i)=C1;
    end
    for i=1:(M3-1)
        A(i,i+1)=D1(i);B(i,i+1)=G1(i);
    end
    for i=1:(M3-2)
        A(i,i+2)=E1(i);B(i,i+2)=-E1(i);
    end
    for i=1:M3
        s(i)=0;
        for j=1:M3
            s(i)=s(i)+B(i,j)*u(j+2,n-2);
        end
        F(i)=s(i);
    end
    Y1=F';
    % Solve a pentadiagonal system Au=F
    d=diag(A);
    a=[diag(A,1);0];
    b=[diag(A,2);0;0];
    c=[0;diag(A,-1)];
    e=[0;0;diag(A,-2)];
    % Solve a pentadiagonal system Au=F
    N=length(Y1);        
    alpha=zeros(N-1,1);beta=zeros(N-2,1);z=zeros(N,1);
    gama=zeros(N,1);mv=zeros(N,1);
    gama(1)=0;     
    % Factor A=LR
    mv(1)=d(1);alpha(1)=a(1)/mv(1);
    beta(1)=b(1)/mv(1);z(1)=Y1(1)/mv(1);
    gama(2)=c(2);
    mv(2)=d(2)-alpha(1)*gama(2);
    alpha(2)=(a(2)-beta(1)*gama(2))/mv(2);
    beta(2)=b(2)/mv(2);
    z(2)=(Y1(2)-z(1)*gama(2))/mv(2);  
    for i=3:N-2
        gama(i)=c(i)-alpha(i-2)*e(i);
        mv(i)=d(i)-beta(i-2)*e(i)-alpha(i-1)*gama(i);
        alpha(i)=(a(i)-beta(i-1)*gama(i))/mv(i);
        beta(i)=b(i)/mv(i);
        z(i)=(Y1(i)-z(i-2)*e(i)-z(i-1)*gama(i))/mv(i);
    end   
    gama(N-1)=c(N-1)-alpha(N-3)*e(N-1);
    mv(N-1)=d(N-1)-beta(N-3)*e(N-1)-alpha(N-2)*gama(N-1);
    alpha(N-1)=(a(N-1)-beta(N-2)*gama(N-1))/mv(N-1);
    gama(N)=c(N)-alpha(N-2)*e(N);
    mv(N)=d(N)-beta(N-2)*e(N)-alpha(N-1)*gama(N);
    z(N-1)=(Y1(N-1)-z(N-3)*e(N-1)-z(N-2)*gama(N-1))/mv(N-1);
    z(N)=(Y1(N)-z(N-2)*e(N)-z(N-1)*gama(N))/mv(N); 
    % Back substitution
    v(N,n)=z(N);v(N-1,n)=z(N-1)-alpha(N-1)*v(N,n);
    for k=N-2:-1:1
        v(k,n)=z(k)-alpha(k)*v(k+1,n)-beta(k)*v(k+2,n);    
    end
    for j=3:(M-1)
        u(j,n)=v(j-2,n);
    end
end
toc;
save u u
S=0;
for j=1:M1
    b(j)=u(j,T1);
    erroru(j)=abs(uexact(j)-b(j));
    S=S+erroru(j)^2;
end
norm2=sqrt(h*S)
norminf=max(erroru)
plot(x,b,'*r',x,uexact,'-k');
ylabel('u');
xlabel('x');
legend('Numerical Result','Exact Solution',2)
axis([-40 60 -0.5 8.5])