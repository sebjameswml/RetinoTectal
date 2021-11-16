% Equations used in my modelling. See tissue.h and
% guidingtissue::exponential_expression and similar

x = [0:0.05:1];

_exp = 0.26 .* exp (2.3.*x) + 1.05;
_exp2 = 0.26 .* exp (1.1.*x) + 1.05;
_log = 2.32 + 1.29 .* log (2.3 .* (x+0.2));
_quad = 1.31 + 2.333 .* x .* x;
_lin = 1.31 + 2.333 .* x;

ki = 1.1
kd = 0.6

_exp_ki = _exp + ki;
_exp_ki_kd = _exp_ki - kd;
_exp_kd = _exp - kd;

_exp2_ki = _exp2 + ki;
_exp2_ki_kd = _exp2_ki - kd;
_exp2_kd = _exp2 - kd;

figure(1); clf;
plot (x, _lin);
hold on;
plot (x, _quad);
plot (x, _log);
plot (x, _exp);
plot (x, _exp2);
lg = legend (['Linear';'Quadratic';'Logarithmic';'Exp';'Exp 2'],'Location','NorthWest');

figure(2); clf;
plot (x, _exp_kd);
hold on;
plot (x, _exp_ki_kd);
plot (x, _exp_ki_kd./_exp_kd);
ylim([0,5]);
lg = legend (['Exp, knocked down';'Exp, knocked down and in';'kikd/kd'],'Location','NorthWest');

figure(3); clf;
plot (x, _exp2_kd);
hold on;
plot (x, _exp2_ki_kd);
plot (x, _exp2_ki_kd./_exp2_kd);
ylim([0,5]);
lg = legend (['Exp2, knocked down';'Exp2, knocked down and in';'kikd/kd'],'Location','NorthWest');

figure(30); clf;
plot (x, _exp2_kd);
hold on;
plot (x, _exp2_ki_kd);
plot (x, _exp2_ki_kd.*_exp2_kd);
ylim([0,5]);
lg = legend (['Exp2, knocked down';'Exp2, knocked down and in';'kikd*kd'],'Location','NorthWest');

% The J effect
figure(4); clf;
plot (x, _exp); % exp for receptors on retina
hold on;
plot (x, flip(_exp)); % ligands
plot (x, _exp .* flip(_exp));
lg = legend (['rcpt';'lgnd';'rcpt*lgnd interaction'],'Location','North');

% with alt exp
figure(14); clf;
plot (x, _exp); % exp for receptors on retina
hold on;
plot (x, _exp2); % exp2 for ligands
plot (x, _exp .* _exp2);
lg = legend (['rcpt';'tec lgnd';'rcpt*lgnd interaction'],'Location','North');

% The I effect (relative)
figure(5); clf;
%_x = -x(2)+x(1);
%x = [_x, x(1:end-1)];
%_exp = 0.26 .* exp (2.3.*x) + 1.05;
plot (x, _exp, '.'); % exp for receptors on retina
hold on;
plot (x, circshift(_exp,1)); % adjacent receptors
plot (x, _exp ./ circshift(_exp,1));
plot (x, circshift(_exp,1) ./ _exp);
lg = legend (['rcpt';'rcptshift';'rcpt/rcptshift interaction';'rcptshift interaction/rcpt'],'Location','North');

% The I effect (mass-action)
figure(6); clf;
plot (x, _exp, 'o-'); % exp for receptors on retina
hold on;
plot (x, circshift(_exp,1), 'o-'); % adjacent receptors
plot (x, _exp .* circshift(_exp,1), 'o-');
lg = legend (['rcpt';'rcptshift';'rcpt*rcptshift interaction'],'Location','North');
