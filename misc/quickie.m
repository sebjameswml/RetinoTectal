A=3.5;
B=0.26;
C=2.3;

A1=0.4+1.05;
B1=0.04;
C1=2.8;
B2=0.048;
C2=3.2;

figure(1)
hold off;
plot (x,  (1.05 + B .* exp(C.*(1-x))));
hold on;

f1 = B1 .* exp (C1.*(1-x));
f2 = B2 .* exp (C2.*(1.2-x));
f = A1+f1+f2;
##plot (x, f);
#plot (x, f1, ':');
#plot (x, f2, ':');
ylim([0,4])

figure(2)
hold off;
plot (x,  (A - (1.05 + B .* exp(C.*(1-x)))));
hold on;
plot (x, (A - f));
ylim([0,4])
