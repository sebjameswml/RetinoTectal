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

% The J effect
figure(1); clf;
plot (x, _exp); % exp for receptors on retina
hold on;
_lig = flip(_exp);
plot (x, _lig); % ligands
plot (x, _exp .* _lig, 'color', colour.orchid2, 'linestyle', '--');
plot (x(1:3), _exp(1) .* _lig(1:3), 'ro-');
plot ([x(1),x(1)], [_exp(2),_exp(1).*_lig(1)], 'ro-');
plot (x(3:7), _exp(5) .* _lig(3:7), 'b');
plot ([x(5),x(5)], [_exp(5),_exp(5).*_lig(5)], 'b-');
plot (x(7:11), _exp(9) .* _lig(7:11), 'g');
plot ([x(9),x(9)], [_exp(9),_exp(9).*_lig(9)], 'g-');
plot (x(11:15), _exp(13) .* _lig(11:15), 'm');
plot ([x(13),x(13)], [_exp(13),_exp(13).*_lig(13)], 'm-');
plot (x(15:19), _exp(17) .* _lig(15:19), 'k');
plot ([x(17),x(17)], [_exp(17),_exp(17).*_lig(17)], 'k-');
lg = legend (['rcpt';'lgnd'],'Location','North');
