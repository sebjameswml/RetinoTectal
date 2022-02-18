% Equations used in my modelling. See tissue.h and
% guidingtissue::exponential_expression and similar

% This models the mass-action receptor-receptor interaction

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

%% The I effect
_r0 = flip(_exp); % One receptor type.
_r2 = _exp;       % A second receptor type.
figure(10); clf;
plot (x, _r0, 'linestyle', '--'); % exp for receptors on retina
hold on;
plot (x, _r2, 'linestyle', '--'); % ligands

plot (x, _r0(1) .* _r0, 'r-');
plot (x, _r0(1) .* _r2, 'r:');

plot (x(3:7), _r0(5) .* _r0(3:7), 'b');
plot (x(3:7), _r0(5) .* _r2(3:7), 'b:');

## plot (x(7:11), _r0(9) .* _r2(7:11), 'g');
## plot ([x(9),x(9)], [_r0(9),_r0(9).*_r2(9)], 'go-');
## plot (x(11:15), _r0(13) .* _r2(11:15), 'm');
## plot ([x(13),x(13)], [_r0(13),_r0(13).*_r2(13)], 'mo-');
## plot (x(15:19), _r0(17) .* _r2(15:19), 'c');
## plot ([x(17),x(17)], [_r0(17),_r0(17).*_r2(17)], 'co-');
## plot (x, _r0(21) .* _r2, 'k');
## plot ([x(21),x(21)], [_r0(21),_r0(21).*_r2(21)], 'ko-');

## Vertical lines
plot ([x(5),x(5)], [_r0(5),_r0(5).*_r0(5)], 'bo-');
plot ([x(1),x(1)], [_r0(1), _r0(1).*_r0(1)], 'ro-');

lg = legend (['r0';'r2';'r0[0] x r0 (self-interaction)';'r0[0] x r2 (related receptor interaction)';'r0[.2] x r0';'r0[.2] x r2'],'Location','North');



strlbl = ['Temporal -------------------------------------> Nasal'];
figure(10)
ylim([0,20])
xlabel(strlbl)
