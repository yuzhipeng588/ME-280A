function output=weakform(N,k)
Gaussian=[
0.000,0.568;
0.538,0.478;
0.906,0.236;
-0.538,0.478;
-0.906,0.237];
E=0.2;
%N=30;
L=1;
k=12;
he=L/N;
J=he/2;
%f=@(x)-1*k^2*sin(2*pi*k*x/L);  
%u=@(x)(-1*L/(4*E*pi^2))*sin(2*pi*k*x/L)+L*x;
%du=@(x)(-1*L*k/(2*E*pi))*cos(2*pi*k*x/L)+L;
%Xkes=@(x,i)J*x+(2*i-1)*L/N;
%GaussianF=@(x)x;
Ke=zeros(3,3,N);
K=zeros(2*N+1,2*N+1);
Re=zeros(3,N);
a=zeros(2*N+1,1);
R=zeros(2*N+1,1);
dthetahat1=@ (x) x-0.5;
dthetahat2=@ (x) -2*x;
dthetahat3=@ (x) x+0.5;
for i=1:N
    Ke(1,:,i)=[integral(@ (x)(x-0.5)*(x-0.5)*E/J,-1,1,'ArrayValued',true),
        integral(@ (x)(x-0.5)*(-2*x)*E/J,-1,1,'ArrayValued',true),
        integral(@ (x)(x-0.5)*(x+0.5)*E/J,-1,1,'ArrayValued',true)];
    Ke(2,2:3,i)=[integral(@ (x)(-2*x)*(-2*x)*E/J,-1,1,'ArrayValued',true),
        integral(@ (x)(-2*x)*(x+0.5)*E/J,-1,1,'ArrayValued',true)];
    Ke(3,3,i)=integral(@ (x)(x+0.5)*(x+0.5)*E/J,-1,1,'ArrayValued',true);
    Ke(2,1,i)=Ke(1,2,i);
    Ke(3,1,i)=Ke(1,3,i);
    Ke(3,2,i)=Ke(2,3,i);
    Re(1,i)=J*Gaussian(1,2)*thetahat1(Gaussian(1,1))*f(Xkes(Gaussian(1,1),i,J,L,N),k,L);
    Re(1,i)=Re(1,i)+J*Gaussian(2,2)*thetahat1(Gaussian(2,1))*f(Xkes(Gaussian(2,1),i,J,L,N),k,L);
    Re(1,i)=Re(1,i)+J*Gaussian(3,2)*thetahat1(Gaussian(3,1))*f(Xkes(Gaussian(3,1),i,J,L,N),k,L);
    Re(1,i)=Re(1,i)+J*Gaussian(4,2)*thetahat1(Gaussian(4,1))*f(Xkes(Gaussian(4,1),i,J,L,N),k,L);
    Re(1,i)=Re(1,i)+J*Gaussian(5,2)*thetahat1(Gaussian(5,1))*f(Xkes(Gaussian(5,1),i,J,L,N),k,L);
    Re(2,i)=J*Gaussian(1,2)*thetahat2(Gaussian(1,1))*f(Xkes(Gaussian(1,1),i,J,L,N),k,L);
    Re(2,i)=Re(2,i)+J*Gaussian(2,2)*thetahat2(Gaussian(2,1))*f(Xkes(Gaussian(2,1),i,J,L,N),k,L);
    Re(2,i)=Re(2,i)+J*Gaussian(3,2)*thetahat2(Gaussian(3,1))*f(Xkes(Gaussian(3,1),i,J,L,N),k,L);
    Re(2,i)=Re(2,i)+J*Gaussian(4,2)*thetahat2(Gaussian(4,1))*f(Xkes(Gaussian(4,1),i,J,L,N),k,L);
    Re(2,i)=Re(2,i)+J*Gaussian(5,2)*thetahat2(Gaussian(5,1))*f(Xkes(Gaussian(5,1),i,J,L,N),k,L);
    Re(3,i)=J*Gaussian(1,2)*thetahat3(Gaussian(1,1))*f(Xkes(Gaussian(1,1),i,J,L,N),k,L);
    Re(3,i)=Re(3,i)+J*Gaussian(2,2)*thetahat3(Gaussian(2,1))*f(Xkes(Gaussian(2,1),i,J,L,N),k,L);
    Re(3,i)=Re(3,i)+J*Gaussian(3,2)*thetahat3(Gaussian(3,1))*f(Xkes(Gaussian(3,1),i,J,L,N),k,L);
    Re(3,i)=Re(3,i)+J*Gaussian(4,2)*thetahat3(Gaussian(4,1))*f(Xkes(Gaussian(4,1),i,J,L,N),k,L);
    Re(3,i)=Re(3,i)+J*Gaussian(5,2)*thetahat3(Gaussian(5,1))*f(Xkes(Gaussian(5,1),i,J,L,N),k,L);
end
for i=2:2*N
    if i==2
        K(1,1)=Ke(1,1,1);
        K(1,2)=Ke(1,2,1);
        K(1,3)=Ke(1,3,1);
        R(1)=Re(1,1);
        K(2,1)=Ke(2,1,1);
        K(2,2)=Ke(2,2,1);
        K(2,3)=Ke(2,3,1);
        R(2)=Re(2,1);
        continue
    end
    if(i==2*N)
        K(2*N,2*N-1)=Ke(2,1,N);
        K(2*N,2*N)=Ke(2,2,N);
        K(2*N,2*N+1)=Ke(2,3,N);
        R(2*N)=Re(2,N);
        K(2*N+1,2*N-1)=Ke(3,1,N);
        K(2*N+1,2*N)=Ke(3,2,N);
        K(2*N+1,2*N+1)=Ke(3,3,N);
        R(2*N+1)=Re(3,N);
        continue
    else
        if(rem (i,2)==1)
            K(i,i-2)=Ke(3,1,floor(i/2));
            K(i,i-1)=Ke(3,2,floor(i/2));
            K(i,i)=Ke(3,3,floor(i/2))+Ke(1,1,floor(i/2)+1);
            K(i,i+1)=Ke(1,2,floor(i/2)+1);
            K(i,i+2)=Ke(1,3,floor(i/2)+1);
            R(i)=Re(3,floor(i/2))+Re(1,floor(i/2)+1);
        end
        if(rem (i,2)==0)
            K(i,i-1)=Ke(2,1,floor(i/2));
            K(i,i)=Ke(2,2,floor(i/2));
            K(i,i+1)=Ke(2,3,floor(i/2));
            R(i)=Re(2,floor(i/2));
        end
    end
end
for i=1:N+1
end

Kc=K(2:2*N,2:2*N);
Rc=R(2:2*N);
Rc(2*N-1)=Rc(2*N-1)+K(2*N,2*N+1);
Rc(2*N-2)=Rc(2*N-2)+K(2*N-1,2*N+1);
Rc(1)=Rc(1)-3*K(2,1);
Rc(2)=Rc(2)-3*K(3,1);
Kc=sparse(Kc);
a=Kc\Rc;
a=[3;a;-1];
x=0:0.002:1;
y1=[];
y2=[];
% figure;
% hold on;
% for i=1:501
%     y1=[y1,u((i-1)*0.002,L,k,E)];
%     if(i==501)
%         y2=[y2,-1];
%         continue;
%     end
%     y2=[y2,uN((i-1)*0.002,he,a,N)];
% end
% plot(x,y1);
% plot(x,y2);
% xlabel('L');
% ylabel('u/uN');
% plot(0:0.05:1,a);
% hold off;
% u=(k*L/(4*pi*pi*E))*(L*sin(2*pi*k*x/L)/(k*pi)-x*cos(2*pi*k*x/L))-(3/pi-4)*x+3;
% du=(k*L/(4*pi*pi*E))*(2*pi*k*x*sin(2*pi*k*x/L)+L*cos(2*pi*k*x/L))-(3/pi-4);
uE=@ (x)E*(((k*L/(4*pi*pi*E))*(2*pi*k*x*sin(2*pi*k*x/L)+L*cos(2*pi*k*x/L))+k/(4*pi*pi*E)-4-((a(2*(floor(x*N)+1)+1)*((2*x-(2*(floor(x*N)+1)-1)*he)/he+0.5)+a(2*(floor(x*N)+1))*(-2*((2*x-(2*(floor(x*N)+1)-1)*he)/he))+a(2*(floor(x*N)+1)-1)*((2*x-(2*(floor(x*N)+1)-1)*he)/he-0.5))/J))^2);
duE=@ (x)E*(((k*L/(4*pi*pi*E))*(2*pi*k*x*sin(2*pi*k*x/L)+L*cos(2*pi*k*x/L))-k/(4*pi*pi*E)-4)^2);
e=(integral(uE,0,L,'ArrayValued',true))^0.5;
uu=(integral(duE,0,L,'ArrayValued',true))^0.5;
eN=e/uu;
output=eN ;
end

