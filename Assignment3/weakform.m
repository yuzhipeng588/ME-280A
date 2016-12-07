%function output=weakform(N)
Gaussian=[
0.000,0.568;
0.538,0.478;
0.906,0.236;
-0.538,0.478;
-0.906,0.237];
N=10000;
L=1;
k=12;
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
for i=1:N
    Ke(:,:,i)=[E1((i-1)*he)/(J*2),-1*E1((i-1)*he)/(J*2);
        -1*E1((i-1)*he)/(J*2),E1((i-1)*he)/(J*2)];
%     cc=(Ke(1,1,i)*Ke(2,2,i))^0.5;
%     Ke(:,:,i)=[1,Ke(1,2,i)/(cc);
%         Ke(2,1,i)/(cc),1];
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
Kc=K(2:N,2:N);
Rc=R(2:N)';
Rc(N-1)=Rc(N-1)-0.7*K(N,N+1);
Rc(1)=Rc(1)+0.3*K(2,1);
T=zeros(N-1,N-1);
KcT=zeros(N-1,N-1);
for i=1:N-1
    T(i,i)=1/sqrt(Kc(i,i));
    KcT(i,i)=1;
    if(i==N-1)
       KcT(i,i-1)=Kc(i,i-1)/(Kc(i,i)*Kc(i-1,i-1))^0.5;
       continue;
    end
    if(i==1)
       KcT(i,i+1)=Kc(i,i+1)/(Kc(i,i)*Kc(i+1,i+1))^0.5;
       continue;
    end
    KcT(i,i-1)=Kc(i,i-1)/(Kc(i,i)*Kc(i-1,i-1))^0.5;
    KcT(i,i+1)=Kc(i,i+1)/(Kc(i,i)*Kc(i+1,i+1))^0.5;
end
KeT=zeros(2,2,N-2);
for i=1:N-2
    KeT(:,:,i)=KcT(i:i+1,i:i+1);
    KeT(2,2,i)=0;
    if(i==N-2)
        KeT(2,2,i)=1;
    end
end
RcT=T'*Rc;
KcT=sparse(KcT);
a=zeros(N-1,1);
j=1;
rT=zeros(N-1,1);
zT=rT;
e=1;
while(e>0.000001)
    apre=a;
    RR=zeros(N-1,1);
    for i=1:N-2
        RR(i:i+1)=RR(i:i+1)+KeT(:,:,i)*a(i:i+1);
    end
    rT=RcT-RR;
    if(j==1)
        zT=rT;
        zKz=0;
        for k=1:N-2
            zKz=zKz+zT(k:k+1)'*KeT(:,:,k)*zT(k:k+1);
        end
        lamda=zT'*rT/zKz;
        a=a+lamda*zT;
    end
    if(j>1)
        zKz=0;
        for k=1:N-2
            zKz=zKz+zT(k:k+1)'*KeT(:,:,k)*zT(k:k+1);
        end
        rKz=0;
        for k=1:N-2
            rKz=rKz+rT(k:k+1)'*KeT(:,:,k)*zT(k:k+1);
        end 
        theta=-rKz/zKz;
        zT=rT+theta*zT;
        zKz=0;
        for k=1:N-2
            zKz=zKz+zT(k:k+1)'*KeT(:,:,k)*zT(k:k+1);
        end
        lamda=zT'*rT/zKz;
        a=a+lamda*zT;
    end
    if(j==1000)
        aeKae=0;
        for k=1:N-2
            aeKae=aeKae+(a(k:k+1)-apre(k:k+1))'*KeT(:,:,k)*(a(k:k+1)-apre(k:k+1));
        end
        ppKpp=0;
        for k=1:N-2
            ppKpp=ppKpp+apre(k:k+1)'*KeT(:,:,k)*apre(k:k+1);
        end        
        e=aeKae/ppKpp;
    end
    j=j+1;
end
% while(e>0.000001)
%     apre=a;
%     for i=1:N-1
%        if(j==1)
%            if(i==1)
%                rT(i:i+1)=RcT(i:i+1)-KcT(i:i+1,i:i+1)*a(i:i+1);
%                zT(i:i+1)=rT(i:i+1);
%                lamda=zT(i:i+1)'*rT(i:i+1)/(zT(i:i+1)'*KcT(i:i+1,i:i+1)*zT(i:i+1));
%                a(i:i+1)=a(i:i+1)+lamda*zT(i:i+1);
% %                a(i)=a(i)+lamda*zT(i);
%                continue;
%            end
%            if(i==N-1)
%                rT(i-1:i)=RcT(i-1:i)-KcT(i,i-1:i)*a(i-1:i);
%                zT(i-1:i)=rT(i-1:i);
%                lamda=zT(i-1:i)'*rT(i-1:i)/(zT(i-1:i)'*KcT(i-1:i,i-1:i)*zT(i-1:i));
%                a(i-1:i)=a(i-1:i)+lamda*zT(i-1:i);
% %                a(i)=a(i)+lamda*zT(i);
%                continue;
%            end
%            rT(i-1:i+1)=RcT(i-1:i+1)-KcT(i-1:i+1,i-1:i+1)*a(i-1:i+1);
%            zT(i-1:i+1)=rT(i-1:i+1);
%            lamda=zT(i-1:i+1)'*rT(i-1:i+1)/(zT(i-1:i+1)'*KcT(i-1:i+1,i-1:i+1)*zT(i-1:i+1));
%            a(i-1:i+1)=a(i-1:i+1)+lamda*zT(i-1:i+1);
% %            a(i)=a(i)+lamda*zT(i);
%        end
%        if(j>1)
%            if(i==1)
%                rT(i:i+1)=RcT(i:i+1)-KcT(i:i+1,i:i+1)*a(i:i+1);
%                the=-rT(i:i+1)'*KcT(i:i+1,i:i+1)*zT(i:i+1)/(zT(i:i+1)'*KcT(i:i+1,i:i+1)*zT(i:i+1));
%                zT(i:i+1)=rT(i:i+1)+the*zT(i:i+1);
%                lamda= zT(i:i+1)'*rT(i:i+1)/(zT(i:i+1)'*KcT(i:i+1,i:i+1)*zT(i:i+1));
% %                if(zT(i:i+1)==[0;0])
% %                    lamda=1;
% %                end
%                a(i:i+1)=a(i:i+1)+lamda*zT(i:i+1);
% %                a(i)=a(i)+lamda*zT(i);
%                continue;
%            end
%            if(i==N-1)
%                rT(i-1:i)=RcT(i-1:i)-KcT(i-1:i,i-1:i)*a(i-1:i);
%                the=-rT(i-1:i)'*KcT(i-1:i,i-1:i)*zT(i-1:i)/(zT(i-1:i)'*KcT(i-1:i,i-1:i)*zT(i-1:i));
%                zT(i-1:i)=rT(i-1:i)+the*zT(i-1:i);
%                lamda= zT(i-1:i)'*rT(i-1:i)/(zT(i-1:i)'*KcT(i-1:i,i-1:i)*zT(i-1:i));
% %                if(zT(i-1:i)==[0;0])
% %                    lamda=1;
% %                end
%                a(i-1:i)=a(i-1:i)+lamda*zT(i-1:i);
% %                a(i)=a(i)+lamda*zT(i);
%                continue;
%            end
%            rT(i-1:i+1)=RcT(i-1:i+1)-KcT(i-1:i+1,i-1:i+1)*a(i-1:i+1);
%            the=-rT(i-1:i+1)'*KcT(i-1:i+1,i-1:i+1)*zT(i-1:i+1)/(zT(i-1:i+1)'*KcT(i-1:i+1,i-1:i+1)*zT(i-1:i+1));
%            zT(i-1:i+1)=rT(i-1:i+1)+the*zT(i-1:i+1);
%            lamda= zT(i-1:i+1)'*rT(i-1:i+1)/(zT(i-1:i+1)'*KcT(i-1:i+1,i-1:i+1)*zT(i-1:i+1));
% %            if(zT(i-1:i+1)==[0;0;0])
% %                lamda=1;
% %            end
%            a(i-1:i+1)=a(i-1:i+1)+lamda*zT(i-1:i+1);
% %            a(i)=a(i)+lamda*zT(i);
%        end
%     end
%     if(mod(j,1000)==0)
%        e=(a-apre)'*KcT*(a-apre)/(apre'*KcT*apre);
%     end
%     j=j+1;
% end
% Kc=sparse(Kc);
% a=Kc\Rc;
% a=[-0.3;a;0.7];
% x=0:0.0002:1;
% y1=[];
% y2=[];
% figure;
% hold on;
% xlabel('L');
% ylabel('u/uN');
% for i=1:5001
% %     y1=[y1,u((i-1)*0.0002,L,k,E)];
%     if(i==5001)
%         y2=[y2,0.7];
%         continue;
%     end
%     y2=[y2,uN((i-1)*0.0002,he,a,N)];
% end
% % plot(x,y1);
% plot(x,y2);

% uE=@(x)E*(((k*L/(4*pi*pi*E))*(2*pi*k*x*sin(2*pi*k*x/L)+L*cos(2*pi*k*x/L))+k/(4*pi*pi*E)-4)-(a(floor(x*N)+2)-a(floor(x*N)+1))/he)^2;
% duE=@(x)E*(((k*L/(4*pi*pi*E))*(2*pi*k*x*sin(2*pi*k*x/L)+L*cos(2*pi*k*x/L))+k/(4*pi*pi*E)-4)^2);
% e=(integral(uE,0,L,'ArrayValued',true))^0.5;
% uu=(integral(duE,0,L,'ArrayValued',true))^0.5;
% eN=e/(integral(duE,0,L,'ArrayValued',true))^0.5;
% output=eN;
% output=a;
%end

% % precondition
% Te=zeros(2,2,N);
% for i=1:N
%     Te(:,:,i)=[1/sqrt(Ke(1,1,i)),0;
%         0,1/sqrt(Ke(2,2,i))];
%     Re(:,i)=Te(:,:,i)*Re(:,i);
% end
% for i=1:N
%     if(i==1)
%         cc=(Ke(1,1,i)*(Ke(2,2,i)+Ke(1,1,i+1)))^0.5;
%     end
%     if(i>1 && i<N)
%         cc=((Ke(1,1,i)+Ke(2,2,i-1))*(Ke(2,2,i)+Ke(1,1,i+1)))^0.5;
%     end
%     if(i==N)
%         cc=((Ke(1,1,i)+Ke(2,2,i-1))*Ke(2,2,i))^0.5;
%     end
%     Ke(:,:,i)=[1,Ke(1,2,i)/(cc);
%     Ke(2,1,i)/(cc),1];
% end
% 
% Ke=Ke(:,:,1:N);
% Re=Re(:,2:N-1);
% Re(2,N-2)=Re(2,N-2)-0.7*Ke(1,2,N);
% Re(1,1)=Re(1,1)+0.3*Ke(2,1,1);
% 
% e=1;
% ae=zeros(2,N-2);
% re=zeros(2,N-2);
% a=zeros(N-1,1);
% ze=re;
% j=1;
% while(e>0.000001)
%     aepre=ae;
%     apre=a;
%     for i=1:N-2
%        if(j==1)
%            re(:,i)=Re(:,i)-Ke(:,:,i+1)*ae(:,i);
%            ze(:,i)=re(:,i);
%            lamda=ze(:,i)'*re(:,i)/(ze(:,i)'*Ke(:,:,i+1)*ze(:,i));
%            ae(:,i)=ae(:,i)+lamda*ze(:,i);
%        end
%        if(j>1)
%            re(:,i)=Re(:,i)-Ke(:,:,i+1)*aepre(:,i);
%            the=-re(:,i)'*Ke(:,:,i+1)*ze(:,i)/(ze(:,i)'*Ke(:,:,i+1)*ze(:,i));
%            ze(:,i)=re(:,i)+the*ze(:,i);
%            lamda= ze(:,i)'*re(:,i)/(ze(:,i)'*Ke(:,:,i+1)*ze(:,i));
%            if(ze(:,i)==[0;0])
%                lamda=1;
%            end
%            ae(:,i)=aepre(:,i)+lamda*ze(:,i);
%        end
%     end
%     j=j+1;
%     e=0;
%     norm = 0;
% %     a=zeros(N-1,1);
% %     a(1)=ae(1,1);
% %     a(N-1)=ae(2,N-2);
% %     K=zeros(N-1,N-1);
% %     K(1,1)=Ke(1,1,2);
% %     K(1,2)=Ke(1,2,2);
% %     K(N-1,N-2)=Ke(2,1,N-1);
% %     K(N-1,N-1)=Ke(2,2,N-1)+Ke(1,1,N);
% %     for i=2:N-2
% %         a(i)=ae(2,i-1)+ae(1,i);
% %         K(i,i-1)=Ke(2,1,i);
% %         K(i,i)=Ke(2,2,i)+Ke(1,1,i+1);
% %         K(i,i+1)=Ke(1,2,i+1);
% %     end
% %     e=(a-apre)'*K*(a-apre)/(apre'*K*apre);
%     for k=1:N-2
%        e=e+(ae(:,k)-aepre(:,k))'*Ke(:,:,k+1)*(ae(:,k)-aepre(:,k));
%        norm = norm + aepre(:,k)'*Ke(:,:,k+1)*aepre(:,k);
%     end
%     e=e/norm;
%end

