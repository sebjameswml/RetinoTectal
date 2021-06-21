%
% This plots the receptor and ligand functions for the fibres
% emanating from the retina in the Gebhardt et al 2012 model; RF and LF.
%
% Also plotted is the ratio of RF and LF and its log, showing that the
% result is a linear function.
%

% Fibre receptors/ligands
alpha = 0.05;
beta = 50;
x = [1:0.01:50];
RF = exp (alpha.*(x - beta./2));
LF = exp (-alpha.*(x - beta./2));

figure(1);
clf;
hold on;
plot (x, RF);
plot (x, LF);
plot (x, LF./RF);
plot (x, log(LF./RF));
plot (x, abs(log(LF./RF)));
legend('RF','LF', 'LF/RF', 'ln(LF/RF)', 'abs(ln(LF/RF))');
title('Retinal fibres')

% Target (tectum) receptor/ligands

% In paper, it's gamma and delta, but I will re-use alpha and beta.
%gamma = 0.05; % same as alpha
%delta = 50; % same as beta

RT = exp (-alpha.*(x - beta./2));
LT = exp (alpha.*(x - beta./2));

figure(2)
clf;
hold on;
plot (x, RT);
plot (x, LT);
plot (x, LT./RT);
plot (x, log(LT./RT));
plot (x, abs(log(LT./RT)));
legend('RT','LT', 'LT/RT', 'ln(LT/RT)', 'abs(ln(LT/RT))');
title('Tectum')
