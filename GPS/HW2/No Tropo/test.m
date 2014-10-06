x=1:5;
y=x;
subplot(2,1,1)
plot(x,y)
hold on
subplot(2,1,2)
y=2*x;
plot(x,y,'r')
legend('First','Second')