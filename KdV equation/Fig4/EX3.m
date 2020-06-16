function EX3
clc;tic;
can1=0.9;can2=1.0e-3;can3=1.0e-5;
XL=0;XR=1;h=0.01;tao=h^2;tfinal=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=ceil((XR-XL)/h);M1=M+1;M5=M-5;
T=ceil(tfinal/tao);T1=T+1;L=XR-XL;
for j=1:M1
    x(j)=XL+(j-1)*h;
end
P2P=5*pi;
for n=1:T1
    t(n)=(n-1)*tao;
    pusai(n)=sin(P2P*t(n));
    dfpusai(n)=P2P*cos(P2P*t(n));
end
% Exact Solution
for j=1:M1
    for n=1:T1
        u(j,n)=0;
        G(j,n)=-(1-x(j)/L)^2*dfpusai(n)+2*can1*pusai(n)*(L-x(j))/(L^2)+2*can2*(pusai(n))^2*(L-x(j))^3/(L^4);
    end
    u(j,1)=-(1-x(j)/L)^2*pusai(1);
    u(j,2)=-(1-x(j)/L)^2*pusai(2);      
end
% Start the iterations
%Boundary conditions
for n=1:T1
    u(1,n)=0;u(2,n)=0;u(3,n)=0;
    u(M-1,n)=0;u(M,n)=0;u(M1,n)=0;
end
for j=1:M5
    for n=1:T1
        v(j,n)=0;
    end
end
%step 2
%solve U2
A=zeros(M5,M5);
A0=can3/(16*h^3);
D0=1/tao;
for j=1:M5
    B0(j)=can1/(24*h)-can3/(2*h^3)+can2/(24*h)*pusai(1)*(1-x(j+3-2)/L)^2;
    C0(j)=-can1/(3*h)+13*can3/(16*h^3)-can2/(3*h)*pusai(1)*(1-x(j+3-1)/L)^2;
    E0(j)=can1/(3*h)-13*can3/(16*h^3)+can2/(3*h)*pusai(1)*(1-x(j+3+1)/L)^2;
    F0(j)=-can1/(24*h)+can3/(2*h^3)-can2/(24*h)*pusai(1)*(1-x(j+3+2)/L)^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=4:M5
    A(j,j-3)=A0;
end
for j=3:M5
    A(j,j-2)=B0(j);
end
for j=2:M5
    A(j,j-1)=C0(j);
end
for j=1:M5
    A(j,j)=D0;
end
for j=1:(M5-1)
    A(j,j+1)=E0(j);
end
for j=1:(M5-2)
    A(j,j+2)=F0(j);
end
for j=1:(M5-3)
    A(j,j+3)=-A0;
end
for j=1:M5
    righ1=-A0*u(j+3-3,1)-B0(j)*u(j+3-2,1)-C0(j)*u(j+3-1,1)+D0*u(j+3,1);
    righ2=-E0(j)*u(j+3+1,1)-F0(j)*u(j+3+2,1)+A0*u(j+3+3,1);
    G0=-can2/(18*h)*(u(j+3+1,2)+u(j+3+1,1)-u(j+3-1,2)-u(j+3-1,1))*(u(j+3+1,2)+u(j+3+1,1)+u(j+3,2)+u(j+3,1)+u(j+3-1,2)+u(j+3-1,1));
    H0=can2/(144*h)*(u(j+3+2,2)+u(j+3+2,1)-u(j+3-2,2)-u(j+3-2,1))*(u(j+3+2,2)+u(j+3+2,1)+u(j+3,2)+u(j+3,1)+u(j+3-2,2)+u(j+3-2,1));
    F(j)=righ1+righ2+G(j+3,1)+G0+H0;
end
N=length(F);  
% Extract bands
c=diag(A);
d=[diag(A,1);0];
e=[diag(A,2);0;0];
f=[diag(A,3);0;0;0];
b=[0;diag(A,-1)];
a=[0;0;diag(A,-2)];
g=[0;0;0;diag(A,-3)];
% Solve a 7-adiagonal system Ax=F
% Extract bands     
alpha=zeros(N,1);beta=zeros(N,1);gama=zeros(N,1);
y=zeros(N,1);
q=zeros(N,1);z=zeros(N,1);mu=zeros(N,1);nv=zeros(N,1);
%j=1
alpha(1)=c(1);beta(1)=d(1)/alpha(1);
q(1)=e(1)/alpha(1);nv(1)=f(1)/alpha(1);
%j=2
gama(2)=b(2);alpha(2)=c(2)-beta(1)*gama(2);
beta(2)=(d(2)-q(1)*gama(2))/alpha(2);
q(2)=(e(2)-nv(1)*gama(2))/alpha(2);
nv(2)=f(2)/alpha(2);
%j=3
z(3)=a(3);gama(3)=b(3)-z(3)*beta(1);
alpha(3)=c(3)-z(3)*q(1)-gama(3)*beta(2);
beta(3)=(d(3)-z(3)*nv(1)-q(2)*gama(3))/alpha(3);
q(3)=(e(3)-nv(2)*gama(3))/alpha(3);
nv(3)=f(3)/alpha(3);
%j=4:N-3
for j=4:N-3
    mu(j)=g(j);
    z(j)=a(j)-mu(j)*beta(j-3);
    gama(j)=b(j)-mu(j)*q(j-3)-z(j)*beta(j-2);
    alpha(j)=c(j)-mu(j)*nv(j-3)-z(j)*q(j-2)-gama(j)*beta(j-1);
    beta(j)=(d(j)-z(j)*nv(j-2)-q(j-1)*gama(j))/alpha(j);
    q(j)=(e(j)-nv(j-1)*gama(j))/alpha(j);
    nv(j)=f(j)/alpha(j);
end
%j=N-2
mu(N-2)=g(N-2);
z(N-2)=a(N-2)-mu(N-2)*beta(N-5);
gama(N-2)=b(N-2)-mu(N-2)*q(N-5)-z(N-2)*beta(N-4);
alpha(N-2)=c(N-2)-mu(N-2)*nv(N-5)-z(N-2)*q(N-4)-gama(N-2)*beta(N-3);
beta(N-2)=(d(N-2)-z(N-2)*nv(N-4)-q(N-3)*gama(N-2))/alpha(N-2);
q(N-2)=(e(N-2)-nv(N-3)*gama(N-2))/alpha(N-2);
%j=N-1
mu(N-1)=g(N-1);
z(N-1)=a(N-1)-mu(N-1)*beta(N-4);
gama(N-1)=b(N-1)-mu(N-1)*q(N-4)-z(N-1)*beta(N-3);
alpha(N-1)=c(N-1)-mu(N-1)*nv(N-4)-z(N-1)*q(N-3)-gama(N-1)*beta(N-2);
beta(N-1)=(d(N-1)-z(N-1)*nv(N-3)-q(N-2)*gama(N-1))/alpha(N-1);
%j=N
mu(N)=g(N);
z(N)=a(N)-mu(N)*beta(N-3);
gama(N)=b(N)-mu(N)*q(N-3)-z(N)*beta(N-2);
alpha(N)=c(N)-mu(N)*nv(N-3)-z(N)*q(N-2)-gama(N)*beta(N-1);
% Back substitution Ly=F
y(1)=F(1)/alpha(1);
y(2)=(F(2)-gama(2)*y(1))/alpha(2);
y(3)=(F(3)-z(3)*y(1)-gama(3)*y(2))/alpha(3);
for j=4:N
    y(j)=(F(j)-mu(j)*y(j-3)-z(j)*y(j-2)-gama(j)*y(j-1))/alpha(j);
end
v(N,2)=y(N);v(N-1,2)=y(N-1)-beta(N-1)*v(N,2);
v(N-2,2)=y(N-2)-beta(N-2)*v(N-1,2)-q(N-2)*v(N,2);
for k=N-3:-1:1
    v(k,2)=y(k)-beta(k)*v(k+1,2)-q(k)*v(k+2,2)-nv(k)*v(k+3,2);    
end
for j=4:(M-2)
    u(j,2)=v(j-3,2);
end
%step 3
Aa=zeros(M5,M5);
A1=can3/(16*h^3);
D1=1/(2*tao);
for n=3:T1
    for j=1:M5
        B1(j)=can1/(24*h)+can2/(72*h)*(u(j+3,n-1)+u(j+3-2,n-1))-can3/(2*h^3)+can2/(24*h)*pusai(n-1)*(1-x(j+3-2)/L)^2;
        C1(j)=-can1/(3*h)-can2/(9*h)*(u(j+3,n-1)+u(j+3-1,n-1))+13*can3/(16*h^3)-can2/(3*h)*pusai(n-1)*(1-x(j+3-1)/L)^2;
        E1(j)=can1/(3*h)+can2/(9*h)*(u(j+3,n-1)+u(j+3+1,n-1))-13*can3/(16*h^3)+can2/(3*h)*pusai(n-1)*(1-x(j+3+1)/L)^2;
        F1(j)=-can1/(24*h)-can2/(72*h)*(u(j+3,n-1)+u(j+3+2,n-1))+can3/(2*h^3)-can2/(24*h)*pusai(n-1)*(1-x(j+3+2)/L)^2;;
    end
    for j=4:M5
        Aa(j,j-3)=A1;
    end
    for j=3:M5
        Aa(j,j-2)=B1(j);
    end
    for j=2:M5
        Aa(j,j-1)=C1(j);
    end
    for j=1:M5
        Aa(j,j)=D1;
    end
    for j=1:(M5-1)
        Aa(j,j+1)=E1(j);
    end
    for j=1:(M5-2)
        Aa(j,j+2)=F1(j);
    end
    for j=1:(M5-3)
        Aa(j,j+3)=-A1;
    end      
    for j=1:M5
        right1=-A1*u(j+3-3,n-2)-B1(j)*u(j+3-2,n-2)-C1(j)*u(j+3-1,n-2)+D1*u(j+3,n-2);
        right2=-E1(j)*u(j+3+1,n-2)-F1(j)*u(j+3+2,n-2)+A1*u(j+3+3,n-2);
        Y(j)=right1+right2+G(j+3,n-1);
    end
    N=length(Y);  
    % Extract bands
    c=diag(Aa);
    d=[diag(Aa,1);0];
    e=[diag(Aa,2);0;0];
    f=[diag(Aa,3);0;0;0];
    b=[0;diag(Aa,-1)];
    a=[0;0;diag(Aa,-2)];
    g=[0;0;0;diag(Aa,-3)];
    % Solve a 7-adiagonal system Ax=Y
    % Extract bands     
    alpha=zeros(N,1);beta=zeros(N,1);gama=zeros(N,1);
    y=zeros(N,1);
    q=zeros(N,1);z=zeros(N,1);mu=zeros(N,1);nv=zeros(N,1);
    %j=1
    alpha(1)=c(1);beta(1)=d(1)/alpha(1);
    q(1)=e(1)/alpha(1);nv(1)=f(1)/alpha(1);
    %j=2
    gama(2)=b(2);alpha(2)=c(2)-beta(1)*gama(2);
    beta(2)=(d(2)-q(1)*gama(2))/alpha(2);
    q(2)=(e(2)-nv(1)*gama(2))/alpha(2);
    nv(2)=f(2)/alpha(2);
    %j=3
    z(3)=a(3);gama(3)=b(3)-z(3)*beta(1);
    alpha(3)=c(3)-z(3)*q(1)-gama(3)*beta(2);
    beta(3)=(d(3)-z(3)*nv(1)-q(2)*gama(3))/alpha(3);
    q(3)=(e(3)-nv(2)*gama(3))/alpha(3);
    nv(3)=f(3)/alpha(3);
    %j=4:N-3
    for j=4:N-3
        mu(j)=g(j);
        z(j)=a(j)-mu(j)*beta(j-3);
        gama(j)=b(j)-mu(j)*q(j-3)-z(j)*beta(j-2);
        alpha(j)=c(j)-mu(j)*nv(j-3)-z(j)*q(j-2)-gama(j)*beta(j-1);
        beta(j)=(d(j)-z(j)*nv(j-2)-q(j-1)*gama(j))/alpha(j);
        q(j)=(e(j)-nv(j-1)*gama(j))/alpha(j);
        nv(j)=f(j)/alpha(j);
    end
    %j=N-2
    mu(N-2)=g(N-2);
    z(N-2)=a(N-2)-mu(N-2)*beta(N-5);
    gama(N-2)=b(N-2)-mu(N-2)*q(N-5)-z(N-2)*beta(N-4);
    alpha(N-2)=c(N-2)-mu(N-2)*nv(N-5)-z(N-2)*q(N-4)-gama(N-2)*beta(N-3);
    beta(N-2)=(d(N-2)-z(N-2)*nv(N-4)-q(N-3)*gama(N-2))/alpha(N-2);
    q(N-2)=(e(N-2)-nv(N-3)*gama(N-2))/alpha(N-2);
    %j=N-1
    mu(N-1)=g(N-1);
    z(N-1)=a(N-1)-mu(N-1)*beta(N-4);
    gama(N-1)=b(N-1)-mu(N-1)*q(N-4)-z(N-1)*beta(N-3);
    alpha(N-1)=c(N-1)-mu(N-1)*nv(N-4)-z(N-1)*q(N-3)-gama(N-1)*beta(N-2);
    beta(N-1)=(d(N-1)-z(N-1)*nv(N-3)-q(N-2)*gama(N-1))/alpha(N-1);
    %j=N
    mu(N)=g(N);
    z(N)=a(N)-mu(N)*beta(N-3);
    gama(N)=b(N)-mu(N)*q(N-3)-z(N)*beta(N-2);
    alpha(N)=c(N)-mu(N)*nv(N-3)-z(N)*q(N-2)-gama(N)*beta(N-1);
    % Back substitution Ly=Y
    y(1)=Y(1)/alpha(1);
    y(2)=(Y(2)-gama(2)*y(1))/alpha(2);
    y(3)=(Y(3)-z(3)*y(1)-gama(3)*y(2))/alpha(3);
    for j=4:N
        y(j)=(Y(j)-mu(j)*y(j-3)-z(j)*y(j-2)-gama(j)*y(j-1))/alpha(j);
    end
    v(N,n)=y(N);v(N-1,n)=y(N-1)-beta(N-1)*v(N,n);
    v(N-2,n)=y(N-2)-beta(N-2)*v(N-1,n)-q(N-2)*v(N,n);
    for k=N-3:-1:1
        v(k,n)=y(k)-beta(k)*v(k+1,n)-q(k)*v(k+2,n)-nv(k)*v(k+3,n);    
    end
    for j=4:(M-2)
        u(j,n)=v(j-3,n);
    end
end
for j=1:M1
    for n=1:T1
        u(j,n)=u(j,n)+(1-x(j)/L)^2*pusai(n);
    end
end
toc;
for n=1:T1
    u1(n)=u(1,n);
    u2(n)=u(M/4+1,n);
    u3(n)=u(M/2+1,n);
    u4(n)=u(M+1,n);
end
plot(t,u1,'-b','LineWidth',2);
xlabel('t');
ylabel('u');
hold on;
plot(t,u2,'-m','LineWidth',2);
hold on;
plot(t,u3,'-r','LineWidth',2);
hold on;
plot(t,u4,'-.k','LineWidth',1);
hold on;
%legend('Forcing','x=L/4','x=L/2','x=L',4)
legend('Forcing','x=L/4','x=L/2',3)
axis([0 1 -1 1])