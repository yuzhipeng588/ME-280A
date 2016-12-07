% All Ke will be recalculated
Gaussian=[0.000000000000000  0.568888888888889;
    0.538469310105683   0.478628670499366;
    0.906179845938664   0.236926885056189;
    -0.538469310105683  0.478628670499366;
    -0.906179845938664  0.236926885056189];   
Ke= zeros(2,2,N);
Re= zeros(2,1,N);
for i= 1:N
    z_x=@(z)nodeIndex(i)*(1-z)/2+nodeIndex(i+1)*(1+z)/2;
    Ke(:,:,i)= E/he(i)*[1,-1;-1,1];
    Re(1,1,i)=Gaussian(:,2)'*(f(z_x(Gaussian(:,1))).*(1-Gaussian(:,1))/2*(he(i)/2));
    Re(2,1,i)=Gaussian(:,2)'*(f(z_x(Gaussian(:,1))).*(1+Gaussian(:,1))/2*(he(i)/2));
end
K=zeros(N+1,N+1);
R=zeros(N+1,1);
for i=1:N
    K(i,i)=K(i,i)+Ke(1,1,i);                   
    K(i,i+1)=K(i,i+1)+Ke(1,2,i);
    K(i+1,i)=K(i+1,i)+Ke(2,1,i);
    K(i+1,i+1)=K(i+1,i+1)+Ke(2,2,i);
    R(i,1)=R(i,1)+Re(1,1,i);
    R(i+1,1)=R(i+1,1)+Re(2,1,i);
end
R(2)=R(2)-u0*K(2,1);
R(N)=R(N)-uL*K(N,N+1);
R=R(2:N,1);
K=K(2:N,2:N);                                       
K=sparse(K);
A=zeros(N+1,1);
A(2:N,1)=K\R;
A(1,1)=u0;
A(N+1,1)=uL;                                        