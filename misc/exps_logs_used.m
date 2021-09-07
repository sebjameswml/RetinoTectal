% Equations used in my modelling. See tissue.h and
% guidingtissue::exponential_expression and similar

x = [0:0.01:1];

_exp = 0.26 .* exp (2.3.*x) + 1.05;
_exp2 = 0.26 .* exp (1.2.*x) + 1.05;
_log = 2.32 + 1.29 .* log (2.3 .* (x+0.2));
_quad = 1.31 + 2.333 .* x .* x;
_lin = 1.31 + 2.333 .* x;

figure(1); clf;
plot (x, _lin);
hold on;
plot (x, _quad);
plot (x, _log);
plot (x, _exp);
plot (x, _exp2);
lg = legend (['Linear';'Quadratic';'Logarithmic';'Exp';'Exp 2'],'Location','NorthWest');
