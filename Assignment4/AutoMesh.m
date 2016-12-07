L=1;
E=1;
u0=1;
uL=1;
f=@(x) 200.*x.^3.*pi.*sin(10.*pi.*x.^5) + 2500.*x.^8.*pi.^2.*cos(10.*pi.*x.^5);
u=@(x)cos(10.*pi.*x.^5);
du=@(x)-50.*x.^4.*pi.*sin(10.*pi.*x.^5);
errorfun1=@(x)E.*du(x).^2;
uE=sqrt(integral(errorfun1,0,1));                        

x=0:0.0002:1;
figure
plot(x,u(x))
hold on
N=20;
nodeIndex=0:L/N:1;              
he=nodeIndex(2:N+1)-nodeIndex(1:N); 
recalculate
AI=[];
for i=1:N
    errorfun=@(x)E*(du(x)-(A(i+1)-A(i))/he(i)).^2;
    AI(i)=sqrt(L/he(i)*integral(errorfun,i*he(i)-he(i),i*he(i)))/uE;        
end
while max(AI)>0.05
    numMesh=find(AI>0.05);       
    nIndex=(nodeIndex(numMesh)+nodeIndex(numMesh+1))/2;
    nodeIndex=[nodeIndex,nIndex];
    nodeIndex=sort(nodeIndex);
    N=size(nodeIndex,2)-1;          
    he=nodeIndex(2:N+1)-nodeIndex(1:N); 
    recalculate
    for i=1:N
        errorfun=@(x)E*(du(x)-(A(i+1)-A(i))/he(i)).^2;
        AI(i)=sqrt(L/he(i)*integral(errorfun,nodeIndex(i),nodeIndex(i+1)))/uE;        
    end      
end
he=[0,he];
theta1=@(x1)0.5*(1-(2*x1-sum(he(1:find(nodeIndex>=x1,1)))-sum(he(1:find(nodeIndex>=x1,1)+1)))/he(find(nodeIndex>=x1,1)+1));
theta2=@(x1)0.5*(1+(2*x1-sum(he(1:find(nodeIndex>=x1,1)))-sum(he(1:find(nodeIndex>=x1,1)+1)))/he(find(nodeIndex>=x1,1)+1));
funuN=@(x1)A(find(nodeIndex>=x1,1))*theta1(x1)+A(find(nodeIndex>=x1,1)+1)*theta2(x1);
uN=[];
for i=1:size(x,2)-1
    uN=[uN,funuN(x(i))];
end
plot(x(1:5000),uN);
num=[];
for i=1:20
    num=[num ,size(find(nodeIndex<0.05*i),2)];
end
num(2:20)=num(2:20)-num(1:19);