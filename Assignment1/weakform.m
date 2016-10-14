function output=weakform(N,k)
Gaussian=[
0.00,0.888;
0.774,0.555;
-0.774,0.555];
E=0.1;
%N=20;
L=1;
%k=1;
he=L/N;
J=he/2;
%f=@(x)-1*k^2*sin(2*pi*k*x/L);  
%u=@(x)(-1*L/(4*E*pi^2))*sin(2*pi*k*x/L)+L*x;
%du=@(x)(-1*L*k/(2*E*pi))*cos(2*pi*k*x/L)+L;
%Xkes=@(x,i)J*x+(2*i-1)*L/N;
%GaussianF=@(x)x;
Ke=zeros(2,2,N);
K=zeros(N+1,N+1);
Re=zeros(2,N);
a=zeros(N+1,1);
R=zeros(N+1,1);
uN=zeros(N,1);
for i=1:N
    Ke(:,:,i)=[E/(J*2),-1*E/(J*2);
        -1*E/(J*2),E/(J*2)];
    Re(1,i)=J*Gaussian(1,2)*thetahat1(Gaussian(1,1))*f(Xkes(Gaussian(1,1),i,J,L,N),k,L);
    Re(1,i)=Re(1,i)+J*Gaussian(2,2)*thetahat1(Gaussian(2,1))*f(Xkes(Gaussian(2,1),i,J,L,N),k,L);
    Re(1,i)=Re(1,i)+J*Gaussian(3,2)*thetahat1(Gaussian(3,1))*f(Xkes(Gaussian(3,1),i,J,L,N),k,L);
    Re(2,i)=J*Gaussian(1,2)*thetahat2(Gaussian(1,1))*f(Xkes(Gaussian(1,1),i,J,L,N),k,L);
    Re(2,i)=Re(2,i)+J*Gaussian(2,2)*thetahat2(Gaussian(2,1))*f(Xkes(Gaussian(2,1),i,J,L,N),k,L);
    Re(2,i)=Re(2,i)+J*Gaussian(3,2)*thetahat2(Gaussian(3,1))*f(Xkes(Gaussian(3,1),i,J,L,N),k,L);
end
for i=1:N+1
    if i==1
        K(i,1)=Ke(1,1,i);
        K(i,2)=Ke(1,2,i);
        R(i)=Re(1,i);
        continue
    end
    if(i==N+1)
        K(i,N)=Ke(2,1,i-1);
        K(i,N+1)=Ke(2,2,i-1);
        R(i)=Re(2,i-1);
        continue
    else
        K(i,i-1)=Ke(2,1,i-1);
        K(i,i)=Ke(2,2,i-1)+Ke(1,1,i);
        K(i,i+1)=Ke(1,2,i);
        R(i)=Re(1,i)+Re(2,i-1);
    end
end
for i=1:N+1
end

Kc=K(2:N,2:N);
Rc=R(2:N);
Rc(N-1)=Rc(N-1)-K(N,N+1);
Kc=sparse(Kc);
a=Kc\Rc;
a=[0;a;1];
x=0:0.01:1;
%figure;
%hold on;
%y1=du(x,L,k,E);
%y2=duN(x,he,a);
%plot(x,y1);
%plot(x,y2);
%plot(0:0.05:1,a);
%hold off;
uE=@(x)E*(((-1*L*k/(2*E*pi))*cos(2*pi*k*x/L)+L)-((a(floor(x*N)+2)-a(floor(x*N)+1)))/he)^2;
duE=@(x)E*(((-1*L*k/(2*E*pi))*cos(2*pi*k*x/L)+L)^2);
e=(integral(uE,0,L,'ArrayValued',true))^0.5;
uu=(integral(duE,0,L,'ArrayValued',true))^0.5;
eN=e/(integral(duE,0,L,'ArrayValued',true))^0.5;
output=eN;
%end

