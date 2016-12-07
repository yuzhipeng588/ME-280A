% NN=[];
%     for n=500:550
%         eN=weakform(n,12);
%         NN=[NN,eN];
%     end    
% a=zeros(N-1,1);
% ae(:,1)=Te(:,:,1)*ae(:,1);
% for i=2:N-2
%     ae(:,i)=Te(:,:,i)*ae(:,i);
%     a(i)=ae(2,i-1)+ae(1,i);
% end
% a(1)=ae(1,1);
% a(N-1)=ae(2,N-2);
% T=zeros(N+1,N+1);
% for i=1:N+1
%     T(i,i)=1/sqrt(K(i,i));
% end
% a=T(2:N,2:N)'*a;
a=T'*a;
a=[-0.3;a;0.7];
x=0:0.00002:1;
y1=[];
y2=[];
figure;
hold on;
xlabel('L');
ylabel('u/uN');
for i=1:50000
%     y1=[y1,du((i-1)*0.0002,L,k,E)];
    y2=[y2,uN((i-1)*0.00002,he,a,N)];
end
% plot(x,y1);
plot(x(1:50000),y2);   
