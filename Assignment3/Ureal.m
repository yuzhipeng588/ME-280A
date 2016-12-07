%Get real solution of 
U =@(x, C, D, EX) D + C.*x - (3.*x.*cos(24.*pi.*x) - sin(24.*pi.*x)./(4.*pi))./(EX.*pi^2);
dU=@(x, C, EX)C + (3*cos(24*pi*x) + 72*x*pi.*sin(24*pi*x))./(EX*pi^2);
syms c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 d1 d2 d3 d4 d5 d6 d7 d8 d9 d10
fun0=U(0,c1,d1,2.5)+0.3;
fun1=U(0.1,c1,d1,2.5)-U(0.1,c2,d2,1);
fun2=U(0.2,c2,d2,1)-U(0.2,c3,d3,1.75);
fun3=U(0.3,c3,d3,1.75)-U(0.3,c4,d4,1.25);
fun4=U(0.4,c4,d4,1.25)-U(0.4,c5,d5,2.75);
fun5=U(0.5,c5,d5,2.75)-U(0.5,c6,d6,3.75);
fun6=U(0.6,c6,d6,3.75)-U(0.6,c7,d7,2.25);
fun7=U(0.7,c7,d7,2.25)-U(0.7,c8,d8,0.75);
fun8=U(0.8,c8,d8,0.75)-U(0.8,c9,d9,2);
fun9=U(0.9,c9,d9,2)-U(0.9,c10,d10,1);
fun10=U(1,c10,d10,1)-0.7;
fun11=2.5.*dU(0.1,c1,2.5)-1.*dU(0.1,c2,1);
fun12=1.*dU(0.2,c2,1)-1.75.*dU(0.2,c3,1.75);
fun13=1.75.*dU(0.3,c3,1.75)-1.25.*dU(0.3,c4,1.25);
fun14=1.25.*dU(0.4,c4,1.25)-2.75.*dU(0.4,c5,2.75);
fun15=2.75.*dU(0.5,c5,2.75)-3.75.*dU(0.5,c6,3.75);
fun16=3.75.*dU(0.6,c6,3.75)-2.25.*dU(0.6,c7,2.25);
fun17=2.25.*dU(0.7,c7,2.25)-0.75.*dU(0.7,c8,0.75);
fun18=0.75.*dU(0.8,c8,0.75)-2.*dU(0.8,c9,2);
fun19=2.*dU(0.9,c9,2)-1.*dU(0.9,c10,1);
c=zeros(1,10);
d=zeros(1,10);
[c(1), c(2), c(3), c(4), c(5), c(6), c(7), c(8), c(9), c(10),...
    d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), d(9), d(10)]=...
    solve(fun0==0,fun1==0,fun2==0,fun3==0,fun4==0,fun5==0,fun6==0,fun7==0,fun8==0,fun9==0,fun10==0,...
    fun11==0,fun12==0,fun13==0,fun14==0,fun15==0,fun16==0,fun17==0,fun18==0,fun19==0,...
    c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10);
E1=[2.5 1.0 1.75 1.25 2.75 3.75 2.25 0.75 2.0 1.0];
E=@(x1)E1(ceil(10.*x1));
u=@(x1)U(x1,c(ceil(10.*x1)+(x1==0)),d(ceil(10.*x1)+(x1==0)),E1(ceil(10.*x1)+(x1==0)));  %real solution of u
du=@(x1)dU(x1,c(ceil(10.*x1)+(x1==0)),E1(ceil(10.*x1)+(x1==0)));
y1=[];
y2=[];
x=0:0.0002:1;
figure;
hold on;
xlabel('L');
ylabel('u/uN');
for i=1:5001
%     y1=[y1,u((i-1)*0.0002,L,k,E)];
    if(i==5001)
        y2=[y2,0.7];
        continue;
    end
    y2=[y2,uN((i-1)*0.0002,he,a,N)];
end
% plot(x,y1);
x1=0:0.001:1;
plot(x,y2);
plot(x1,u(x1));
legend('uN ','uReal');
hold off;
uE=@(x)E1(ceil(10.*x)).*((du(x))-(a(floor(x*N)+2)-a(floor(x*N)+1))/he)^2;
duE=@(x)E1(ceil(10.*x))*(du(x))^2;
e=(integral(uE,0,L,'ArrayValued',true))^0.5;
eN1=e/(integral(duE,0,L,'ArrayValued',true))^0.5;
k=12;
L=1;
uu=@(x)(a(floor(x*N)+2)-a(floor(x*N)+1))/he*E1(ceil(10.*x))*(a(floor(x*N)+2)-a(floor(x*N)+1))/he;
fx=@(x)(-x*(k^3)*cos(2*pi*k*x/L))*(a(floor(x*N)+1)*theta1(x,he,N)+a(floor(x*N)+1+1)*theta2(x,he,N));
J=0.5*integral(uu,0,L,'ArrayValued',true)-integral(fx,0,L,'ArrayValued',true);
% N=1000;
% a=weakform(N);
% uE=@(x)E1(ceil(10.*x)).*((du(x))-(a(floor(x*N)+2)-a(floor(x*N)+1))/he)^2;
% duE=@(x)E1(ceil(10.*x))*(du(x))^2;
% e=(integral(uE,0,L,'ArrayValued',true))^0.5;
% eN1=e/(integral(duE,0,L,'ArrayValued',true))^0.5;
% eN=[eN,eN1];
% N=10000;
% a=weakform(N);
% uE=@(x)E1(ceil(10.*x)).*((du(x))-(a(floor(x*N)+2)-a(floor(x*N)+1))/he)^2;
% duE=@(x)E1(ceil(10.*x))*(du(x))^2;
% e=(integral(uE,0,L,'ArrayValued',true))^0.5;
% eN1=e/(integral(duE,0,L,'ArrayValued',true))^0.5;
% eN=[eN,eN1];
% U={@(x)u(x,c(1),d(1),2.5),@(x)u(x,c(2),d(2),1),@(x)u(x,c(3),d(3),1.75),@(x)u(x,c(4),d(4),1.25),@(x)u(x,c(5),d(5),2.75)...
%     @(x)u(x,c(6),d(6),3.75),@(x)u(x,c(7),d(7),2.25),@(x)u(x,c(8),d(8),0.75),@(x)u(x,c(9),d(9),2),@(x)u(x,c(10),d(10),1)};
% hold on
% x1=0.1:0.001:0.2;
% plot(x1,U{2}(x1))
% x1=0.2:0.001:0.3;
% plot(x1,U{3}(x1))
% x1=0.3:0.001:0.4;
% plot(x1,U{4}(x1))
% x1=0.4:0.001:0.5;
% plot(x1,U{5}(x1))
% x1=0.5:0.001:0.6;
% plot(x1,U{6}(x1))
% x1=0.6:0.001:0.7;
% plot(x1,U{7}(x1))
% x1=0.7:0.001:0.8;
% plot(x1,U{8}(x1))
% x1=0.8:0.001:0.9;
% plot(x1,U{9}(x1))
% x1=0.9:0.001:1;
% plot(x1,U{10}(x1))
% 
% 
% 
% % if
% %     x>=0 & x<=0.1
% %     U=@(x)u(x,c(1),d(1),2.5);
% % elseif x>0.1 & x<=0.2
% %     U=@(x)=u(x,c(2),d(2),1)+...
% % elseif x>0.2 & x<=0.3
% %     U=@(x)u(x,c(3),d(3),1.75)+...
% % elseif x>0.3 & x<=0.4
% %     U=@(x)u(x,c(4),d(4),1.25)+...
% %     (x>0.4 & x<=0.5)..*u(x,c(5),d(5),2.75)+...
% %     (x>0.5 & x<=0.6)..*u(x,c(6),d(6),3.75)+...
% %     (x>0.6 & x<=0.7)..*u(x,c(7),d(7),2.25)+...
% %     (x>0.7 & x<=0.8)..*u(x,c(8),d(8),0.75)+...
% %     (x>0.8 & x<=0.9)..*u(x,c(9),d(9),2)+...
% %     (x>0.9 & x<=1)..*u(x,c(10),d(10),1);
% % end