function EX3ok
clc;tic;
can1=1;can2=1;can3=1;
XL=0;XR=1;h=0.001;tao=0.0005;tfinal=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=ceil((XR-XL)/h);M1=M+1;M5=M-5;
T=ceil(tfinal/tao);T1=T+1;L=XR-XL;
for j=1:M1
    x(j)=XL+(j-1)*h;
end
for n=1:T1
    t(n)=(n-1)*tao;
    pusai(n)=0.1*cos(2*pi*t(n));
    dfpusai(n)=-0.2*pi*sin(2*pi*t(n));
    dfdfpusai(n)=-0.2*pi*2*pi*cos(2*pi*t(n));
end
% Exact Solution
for j=1:M1
    for n=1:T1
        u(j,n)=0;
        G(j,n)=-(1-x(j)/L)^2*dfpusai(n)+2*can1*pusai(n)*(L-x(j))/(L^2)+2*can2*((pusai(n))^2)*((L-x(j))^3)/(L^4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        uexa00=pusai(n)*(1-x(j))^2;
        uexa11=pusai(n)/12*((1-x(j))^2-(1-x(j))^4);
        uexa12=(((pusai(n))^3)/30-dfpusai(n)/60)*((1-x(j))^2-(1-x(j))^5);
        uexa21=1/24*(pusai(n)/6+((pusai(n))^3)/15-dfpusai(n)/30)*((1-x(j))^2-(1-x(j))^4);
        uexa22=1/24*(dfdfpusai(n)/60-((pusai(n))^2)*dfpusai(n)/10-dfpusai(n)/12)*((1-x(j))^2-(1-x(j))^5);
        uexa23=1/120*(((pusai(n))^2)/3-pusai(n)/3-pusai(n)*dfpusai(n)/15+2*((pusai(n))^4)/15)*((1-x(j))^2-(1-x(j))^6);
        uexa24=1/210*(pusai(n)/6-((pusai(n))^3)/6)*((1-x(j))^2-(1-x(j))^7);
        uexa25=1/336*(-((pusai(n))^2)/2-dfdfpusai(n)/60+dfpusai(n)*((pusai(n))^2)/10)*((1-x(j))^2-(1-x(j))^8);
        uexa26=1/504*(-7*((pusai(n))^4)/30+7*pusai(n)*dfpusai(n)/60)*((1-x(j))^2-(1-x(j))^9);
        uexat(j,n)=uexa00+uexa11+uexa12+uexa21+uexa22+uexa23+uexa24+uexa25+uexa26;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    ADM(j)=uexat(j,T1);
    erroru(j)=0;
end
% Start the iterations
%Boundary conditions
for n=1:T1
    u(1,n)=0;u(2,n)=0;u(3,n)=0;
    u(M-1,n)=0;u(M,n)=0;u(M1,n)=0;
    for j=1:M5
        v(j,n)=0;
    end
end
%step 2  %solve U2
A=zeros(M5,M5);A0=can3/(16*h^3);D0=1/tao;
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
    righ=-A0*u(j+3-3,1)-B0(j)*u(j+3-2,1)-C0(j)*u(j+3-1,1)+D0*u(j+3,1)-E0(j)*u(j+3+1,1)-F0(j)*u(j+3+2,1)+A0*u(j+3+3,1);
    G0=-can2/(18*h)*(u(j+3+1,2)+u(j+3+1,1)-u(j+3-1,2)-u(j+3-1,1))*(u(j+3+1,2)+u(j+3+1,1)+u(j+3,2)+u(j+3,1)+u(j+3-1,2)+u(j+3-1,1));
    H0=can2/(144*h)*(u(j+3+2,2)+u(j+3+2,1)-u(j+3-2,2)-u(j+3-2,1))*(u(j+3+2,2)+u(j+3+2,1)+u(j+3,2)+u(j+3,1)+u(j+3-2,2)+u(j+3-2,1));
    F0(j)=righ+G(j+3,1)+G0+H0;
end
N=length(F0);  
% Extract bands
c=diag(A);
d=[diag(A,1);0];
e=[diag(A,2);0;0];
f=[diag(A,3);0;0;0];
b=[0;diag(A,-1)];
a=[0;0;diag(A,-2)];
g=[0;0;0;diag(A,-3)];
% Solve a 7-adiagonal system Ax=F0
% Extract bands     
alpha=zeros(N,1);beta=zeros(N,1);gama=zeros(N,1);y=zeros(N,1);
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
y(1)=F0(1)/alpha(1);
y(2)=(F0(2)-gama(2)*y(1))/alpha(2);
y(3)=(F0(3)-z(3)*y(1)-gama(3)*y(2))/alpha(3);
for j=4:N
    y(j)=(F0(j)-mu(j)*y(j-3)-z(j)*y(j-2)-gama(j)*y(j-1))/alpha(j);
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
B=zeros(M5,M5);A1=can3/(16*h^3);D1=1/(2*tao);
for n=3:T1
    for j=1:M5
        B1(j)=can1/(24*h)+can2/(72*h)*(u(j+3,n-1)+u(j+3-2,n-1))-can3/(2*h^3)+can2/(24*h)*pusai(n-1)*((1-x(j+3-2)/L)^2);
        C1(j)=-can1/(3*h)-can2/(9*h)*(u(j+3,n-1)+u(j+3-1,n-1))+13*can3/(16*h^3)-can2/(3*h)*pusai(n-1)*((1-x(j+3-1)/L)^2);
        E1(j)=can1/(3*h)+can2/(9*h)*(u(j+3,n-1)+u(j+3+1,n-1))-13*can3/(16*h^3)+can2/(3*h)*pusai(n-1)*((1-x(j+3+1)/L)^2);
        F1(j)=-can1/(24*h)-can2/(72*h)*(u(j+3,n-1)+u(j+3+2,n-1))+can3/(2*h^3)-can2/(24*h)*pusai(n-1)*((1-x(j+3+2)/L)^2);
    end
    for j=4:M5
        B(j,j-3)=A1;
    end
    for j=3:M5
        B(j,j-2)=B1(j);
    end
    for j=2:M5
        B(j,j-1)=C1(j);
    end
    for j=1:M5
        B(j,j)=D1;
    end
    for j=1:(M5-1)
        B(j,j+1)=E1(j);
    end
    for j=1:(M5-2)
        B(j,j+2)=F1(j);
    end
    for j=1:(M5-3)
        B(j,j+3)=-A1;
    end      
    for j=1:M5
        riht=-A1*u(j+3-3,n-2)-B1(j)*u(j+3-2,n-2)-C1(j)*u(j+3-1,n-2)+D1*u(j+3,n-2)-E1(j)*u(j+3+1,n-2)-F1(j)*u(j+3+2,n-2)+A1*u(j+3+3,n-2);
        Y(j)=riht+G(j+3,n-1);
    end
    N=length(Y);  
    % Extract bands
    c=diag(B);
    d=[diag(B,1);0];
    e=[diag(B,2);0;0];
    f=[diag(B,3);0;0;0];
    b=[0;diag(B,-1)];
    a=[0;0;diag(B,-2)];
    g=[0;0;0;diag(B,-3)];
    % Solve a 7-adiagonal system Ax=Y
    % Extract bands     
    alpha=zeros(N,1);beta=zeros(N,1);gama=zeros(N,1);
    y=zeros(N,1);q=zeros(N,1);z=zeros(N,1);mu=zeros(N,1);nv=zeros(N,1);
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    % Back substitution Lyy=Y
    yy(1)=Y(1)/alpha(1);
    yy(2)=(Y(2)-gama(2)*yy(1))/alpha(2);
    yy(3)=(Y(3)-z(3)*yy(1)-gama(3)*yy(2))/alpha(3);
    for j=4:N
        yy(j)=(Y(j)-mu(j)*yy(j-3)-z(j)*yy(j-2)-gama(j)*yy(j-1))/alpha(j);
    end
    v(N,n)=yy(N);
    v(N-1,n)=yy(N-1)-beta(N-1)*v(N,n);
    v(N-2,n)=yy(N-2)-beta(N-2)*v(N-1,n)-q(N-2)*v(N,n);
    for k=N-3:-1:1
        v(k,n)=yy(k)-beta(k)*v(k+1,n)-q(k)*v(k+2,n)-nv(k)*v(k+3,n);    
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
save u u
toc;
S=0;
for j=1:M1
    b(j)=u(j,T1);
    erroru(j)=abs(ADM(j)-b(j));
    S=S+erroru(j)^2;
end
norm2=sqrt(h*S)
norminf=max(erroru)
plot(x,b,'.r',x,ADM,'-k');
ylabel('u');
xlabel('x');
%axis([-40 100 -0.01 0.51])
legend('Present scheme','ADM method [31]',2)