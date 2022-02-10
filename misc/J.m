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

%% The J effect
_r0 = _exp;
_l0 = flip(_exp);
figure(1); clf;
plot (x, _r0); % exp for receptors on retina
hold on;
plot (x, _l0); % ligands
plot (x, _r0 .* _l0, 'color', colour.orchid2, 'linestyle', '--');
plot (x(1:3), _r0(1) .* _l0(1:3), 'ro-');
plot ([x(1),x(1)], [_r0(2),_r0(1).*_l0(1)], 'ro-');
plot (x(3:7), _r0(5) .* _l0(3:7), 'b');
plot ([x(5),x(5)], [_r0(5),_r0(5).*_l0(5)], 'b-');
plot (x(7:11), _r0(9) .* _l0(7:11), 'g');
plot ([x(9),x(9)], [_r0(9),_r0(9).*_l0(9)], 'g-');
plot (x(11:15), _r0(13) .* _l0(11:15), 'm');
plot ([x(13),x(13)], [_r0(13),_r0(13).*_l0(13)], 'm-');
plot (x(15:19), _r0(17) .* _l0(15:19), 'k');
plot ([x(17),x(17)], [_r0(17),_r0(17).*_l0(17)], 'k-');
lg = legend (['rcpt0';'lgnd0'],'Location','North');

_r1 = flip(_exp); % opposing receptors
_l1 = _exp; % opposing ligands
figure(2); clf;
plot (x, _r1);
hold on;
plot (x, _l1); % ligands
plot (x, _r1 .* _l1, 'color', colour.orchid2, 'linestyle', '--');
plot (x(1:3), _r1(1) .* _l1(1:3), 'ro-');
plot ([x(1),x(1)], [_r1(2),_r1(1).*_l1(1)], 'ro-');
plot (x(3:7), _r1(5) .* _l1(3:7), 'b');
plot ([x(5),x(5)], [_r1(5),_r1(5).*_l1(5)], 'b-');
plot (x(7:11), _r1(9) .* _l1(7:11), 'g');
plot ([x(9),x(9)], [_r1(9),_r1(9).*_l1(9)], 'g-');
plot (x(11:15), _r1(13) .* _l1(11:15), 'm');
plot ([x(13),x(13)], [_r1(13),_r1(13).*_l1(13)], 'm-');
plot (x(15:19), _r1(17) .* _l1(15:19), 'k');
plot ([x(17),x(17)], [_r1(17),_r1(17).*_l1(17)], 'k-');
lg = legend (['rcpt1';'lgnd1'],'Location','North');

%% Now plot combined repulsion due to both receptor/ligand pairs
figure(1)
plot (x, _r1(1) .* _l1 + _r0(1) .* _l0, 'r--');
plot (x, _r1(5) .* _l1 + _r0(5) .* _l0, 'b--');
plot (x, _r1(9) .* _l1 + _r0(9) .* _l0, 'g--');
plot (x, _r1(13) .* _l1 + _r0(13) .* _l0, 'm--');
plot (x, _r1(17) .* _l1 + _r0(17) .* _l0, 'k--');
ylim([0,20])
figure(2)
plot (x(1:3), _r1(1) .* _l1(1:3) + _r0(1) .* _l0(1:3), 'r--');
plot (x(3:7), _r1(5) .* _l1(3:7) + _r0(5) .* _l0(3:7), 'b--');
plot (x(7:11), _r1(9) .* _l1(7:11) + _r0(9) .* _l0(7:11), 'g--');
plot (x(11:15), _r1(13) .* _l1(11:15) + _r0(13) .* _l0(11:15), 'm--');
plot (x(15:19), _r1(17) .* _l1(15:19) + _r0(17) .* _l0(15:19), 'k--');
ylim([0,20])
