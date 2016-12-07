function output=weakform(N,k)
Gaussian=[
0.000,0.568;
0.538,0.478;
0.906,0.236;
-0.538,0.478;
-0.906,0.237];
N=48;
E=0.2;
L=1;
k=12;
he=L/N;
J=he/2;
%f=@(x)-1*k^2*sin(2*pi*k*x/L);  
%u=@(x)(-1*L/(4*E*pi^2))*sin(2*pi*k*x/L)+L*x;
%du=@(x)(-1*L*k/(2*E*pi))*cos(2*pi*k*x/L)+L;
%Xkes=@(x,i)J*x+(2*i-1)*L/N;
%GaussianF=@(x)x;
Ke=zeros(4,4,N);
K=zeros(3*N+1,3*N+1);
Re=zeros(4,N);
a=zeros(3*N+1,1);
R=zeros(3*N+1,1);
% thetahat1=@ (x) (-9/16)*(x+1/3)*(x-1/3)*(x-1);
dthetahat1=@ (x) (-9/16)*(3*x*x-2*x-1/9);
dthetahat2=@ (x) (27/16)*(3*x*x-x*2/3-1);
dthetahat3=@ (x) (-27/16)*(3*x*x+x*2/3-1);
dthetahat4=@ (x) (9/16)*(3*x*x+2*x-1/9);
for i=1:N
    Ke(1,:,i)=[integral(@ (x)(-9/16)*(3*x*x-2*x-1/9)*(-9/16)*(3*x*x-2*x-1/9)*E/J,-1,1,'ArrayValued',true),
        integral(@ (x)(-9/16)*(3*x*x-2*x-1/9)*(27/16)*(3*x*x-x*2/3-1)*E/J,-1,1,'ArrayValued',true),
        integral(@ (x)(-9/16)*(3*x*x-2*x-1/9)*(-27/16)*(3*x*x+x*2/3-1)*E/J,-1,1,'ArrayValued',true),
        integral(@ (x)(-9/16)*(3*x*x-2*x-1/9)*(9/16)*(3*x*x+2*x-1/9)*E/J,-1,1,'ArrayValued',true)];
    Ke(2,2:4,i)=[integral(@ (x)(27/16)*(3*x*x-x*2/3-1)*(27/16)*(3*x*x-x*2/3-1)*E/J,-1,1,'ArrayValued',true),
        integral(@ (x)(27/16)*(3*x*x-x*2/3-1)*(-27/16)*(3*x*x+x*2/3-1)*E/J,-1,1,'ArrayValued',true),
        integral(@ (x)(27/16)*(3*x*x-x*2/3-1)*(9/16)*(3*x*x+2*x-1/9)*E/J,-1,1,'ArrayValued',true)];
    Ke(3,3:4,i)=[integral(@ (x)(-27/16)*(3*x*x+x*2/3-1)*(-27/16)*(3*x*x+x*2/3-1)*E/J,-1,1,'ArrayValued',true),
        integral(@ (x)(-27/16)*(3*x*x+x*2/3-1)*(9/16)*(3*x*x+2*x-1/9)*E/J,-1,1,'ArrayValued',true)];
    Ke(4,4,i)=integral(@ (x)(9/16)*(3*x*x+2*x-1/9)*(9/16)*(3*x*x+2*x-1/9)*E/J,-1,1,'ArrayValued',true);
    Ke(2,1,i)=Ke(1,2,i);
    Ke(3,1,i)=Ke(1,3,i);
    Ke(4,1,i)=Ke(1,4,i);
    Ke(3,2,i)=Ke(2,3,i);
    Ke(4,2,i)=Ke(2,4,i);
    Ke(4,3,i)=Ke(3,4,i);
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
    Re(4,i)=J*Gaussian(1,2)*thetahat4(Gaussian(1,1))*f(Xkes(Gaussian(1,1),i,J,L,N),k,L);
    Re(4,i)=Re(4,i)+J*Gaussian(2,2)*thetahat4(Gaussian(2,1))*f(Xkes(Gaussian(2,1),i,J,L,N),k,L);
    Re(4,i)=Re(4,i)+J*Gaussian(3,2)*thetahat4(Gaussian(3,1))*f(Xkes(Gaussian(3,1),i,J,L,N),k,L);
    Re(4,i)=Re(4,i)+J*Gaussian(4,2)*thetahat4(Gaussian(4,1))*f(Xkes(Gaussian(4,1),i,J,L,N),k,L);
    Re(4,i)=Re(4,i)+J*Gaussian(5,2)*thetahat4(Gaussian(5,1))*f(Xkes(Gaussian(5,1),i,J,L,N),k,L);
end
for i=3:3*N-1
    if i==3
        K(1,1)=Ke(1,1,1);
        K(1,2)=Ke(1,2,1);
        K(1,3)=Ke(1,3,1);
        K(1,4)=Ke(1,4,1);
        R(1)=Re(1,1);
        K(2,1)=Ke(2,1,1);
        K(2,2)=Ke(2,2,1);
        K(2,3)=Ke(2,3,1);
        K(2,4)=Ke(2,4,1);
        R(2)=Re(2,1);
        K(3,1)=Ke(3,1,1);
        K(3,2)=Ke(3,2,1);
        K(3,3)=Ke(3,3,1);
        K(3,4)=Ke(3,4,1);
        R(3)=Re(3,1);
        continue
    end
    if(i==3*N-1)
        K(3*N-1,3*N-2)=Ke(2,1,N);
        K(3*N-1,3*N-1)=Ke(2,2,N);
        K(3*N-1,3*N)=Ke(2,3,N);
        K(3*N-1,3*N+1)=Ke(2,4,N);
        R(3*N-1)=Re(2,N);
        K(3*N,3*N-2)=Ke(3,1,N);
        K(3*N,3*N-1)=Ke(3,2,N);
        K(3*N,3*N)=Ke(3,3,N);
        K(3*N,3*N+1)=Ke(3,4,N);
        R(3*N)=Re(3,N);
        K(3*N+1,3*N-2)=Ke(4,1,N);
        K(3*N+1,3*N-1)=Ke(4,2,N);
        K(3*N+1,3*N)=Ke(4,3,N);
        K(3*N+1,3*N+1)=Ke(4,4,N);
        R(3*N+1)=Re(4,N);
        continue
    else
        if(rem (i,3)==1)
            K(i,i-3)=Ke(4,1,floor(i/3));
            K(i,i-2)=Ke(4,2,floor(i/3));
            K(i,i-1)=Ke(4,3,floor(i/3));
            K(i,i)=Ke(4,4,floor(i/3))+Ke(1,1,floor(i/3)+1);
            K(i,i+1)=Ke(1,2,floor(i/3)+1);
            K(i,i+2)=Ke(1,3,floor(i/3)+1);
            K(i,i+3)=Ke(1,4,floor(i/3)+1);
            R(i)=Re(4,floor(i/3))+Re(1,floor(i/3)+1);
        end
        if(rem (i,3)==2)
            K(i,i-1)=Ke(2,1,floor(i/3));
            K(i,i)=Ke(2,2,floor(i/3));
            K(i,i+1)=Ke(2,3,floor(i/3));
            K(i,i+2)=Ke(2,4,floor(i/3));
            R(i)=Re(2,floor(i/3)+1);
        end
        if(rem (i,3)==0)
            K(i,i-2)=Ke(3,1,floor(i/3));
            K(i,i-1)=Ke(3,2,floor(i/3));
            K(i,i)=Ke(3,3,floor(i/3));
            K(i,i+1)=Ke(3,4,floor(i/3));
            R(i)=Re(3,floor(i/3));
        end
    end
end
for i=1:N+1
end

Kc=K(2:3*N,2:3*N);
Rc=R(2:3*N);
Rc(3*N-1)=Rc(3*N-1)+K(3*N,3*N+1);
Rc(3*N-2)=Rc(3*N-2)+K(3*N-1,3*N+1);
Rc(3*N-3)=Rc(3*N-3)+K(3*N-2,3*N+1);
Rc(1)=Rc(1)-3*K(2,1);
Rc(2)=Rc(2)-3*K(3,1);
Rc(3)=Rc(3)-3*K(4,1);
Kc=sparse(Kc);
a=Kc\Rc;
a=[3;a;-1];
x=0:0.002:1;
y1=[];
y2=[];
figure;
hold on;
for i=1:501
    y1=[y1,u((i-1)*0.002,L,k,E)];
    if(i==501)
        y2=[y2,-1];
        continue;
    end
    y2=[y2,uN((i-1)*0.002,he,a,N)];
end
plot(x,y1);
plot(x,y2);
xlabel('L');
ylabel('u/uN');
hold off;
%i=(floor(x*N)+1);
uE=@ (x)E*(((k*L/(4*pi*pi*E))*(2*pi*k*x*sin(2*pi*k*x/L)+L*cos(2*pi*k*x/L))+k/(4*pi*pi*E)-4-(a(3*(floor(x*N)+1)-2)*((-9/16)*(3*((2*x-(2*(floor(x*N)+1)-1)*he)/he)*((2*x-(2*(floor(x*N)+1)-1)*he)/he)-2*((2*x-(2*(floor(x*N)+1)-1)*he)/he)-1/9))+a(3*(floor(x*N)-1))*((27/16)*(3*((2*x-(2*(floor(x*N)+1)-1)*he)/he)*((2*x-(2*(floor(x*N)+1)-1)*he)/he)-((2*x-(2*(floor(x*N)+1)-1)*he)/he)*2/3-1))+a(3*(floor(x*N)+1))*((-27/16)*(3*((2*x-(2*(floor(x*N)+1)-1)*he)/he)*((2*x-(2*(floor(x*N)+1)-1)*he)/he)+((2*x-(2*(floor(x*N)+1)-1)*he)/he)*2/3-1))+a(3*(floor(x*N)+1)+1)*((9/16)*(3*((2*x-(2*(floor(x*N)+1)-1)*he)/he)*((2*x-(2*(floor(x*N)+1)-1)*he)/he)+2*((2*x-(2*(floor(x*N)+1)-1)*he)/he)-1/9)))/J)^2);
duE=@ (x)E*(((k*L/(4*pi*pi*E))*(2*pi*k*x*sin(2*pi*k*x/L)+L*cos(2*pi*k*x/L))-k/(4*pi*pi*E)-4)^2);
e=(integral(uE,0,L,'ArrayValued',true))^0.5;
uu=(integral(duE,0,L,'ArrayValued',true))^0.5;
eN=e/uu;
output=eN ;
end

