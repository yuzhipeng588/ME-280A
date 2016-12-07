f1=figure(1);
set(f1,'name','N for each k','Numbertitle','off');
subplot(4,2,1);
hold on;
plot([1,2,4,8,16,32], eN(:,5));
%ylim([0,0.1]);
legend('error');
xlabel('k');
ylabel('error');

subplot(4,1,3);
plot(2:4096, eN(3,2:4096));
%ylim([0,0.1]);
%xlim([2,30]);
legend('error');
xlabel('N');
ylabel('error');

subplot(4,1,2);
plot(2:4096, eN(2,2:4096));
%ylim([0,0.1]);
%xlim([2,30]);
legend('error');
xlabel('N');
ylabel('error');

subplot(4,1,4);
plot(2:4096, eN(1,2:4096));
%ylim([0,0.1]);
%xlim([2,30]);
legend('error');
xlabel('N');
ylabel('error');
box on;
hold off;