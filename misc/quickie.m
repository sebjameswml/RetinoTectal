x=[0:0.01:1];
epha4_constant=3.5;
B=0.26;
C=2.3;
D=1.05;

ephrinA = D + B .* exp(C.*(1-x));

w_EphAx = 0.3;
w_EphA4 = w_EphAx .* 0.611;

A1=0.83;
B1=0.04;
C1=2.8;
B2=0.048;
C2=2.8;

figure(1)
hold off;
_p_epha4 = w_EphA4 .* epha4_constant * ephrinA;
plot (x, _p_epha4, 'g-');
hold on;

f1 = B1 .* exp (C1.*(1-x));
f2 = B2 .* exp (C2.*(1.2-x));
f = A1+f1+f2;
plot (x, f, 'c-');
plot (x, f1, ':');
plot (x, f2, ':');
ylim([0,4])

figure(2)
hold off;
plot (x,  (epha4_constant - _p_epha4), 'g-');
hold on;
plot (x, (epha4_constant - f), 'c-');

A3=0.95;
B3=0.04;
C3=2.8;
B4=0.048;
C4=3;

f3 = B3 .* exp (C3.*(1-x));
f4 = B4 .* exp (C4.*(1.2-x));
f5 = A3+f3+f4;
plot (x, (epha4_constant - f5), 'r-');
legend(['WT';'kd1';'kd2'])

ylim([0,4])
