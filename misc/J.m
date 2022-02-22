% Equations used in my modelling. See tissue.h and
% guidingtissue::exponential_expression and similar

% Set up the one figure
h_f = figure(1); clf;
h_f_pos = get(h_f, 'Position');
w=2500;
h=1280;
set(h_f, 'Position', [20, 1500, w, h]);
clf;

## What ligand expression profile?
ligand_expression = 'lin'; #  or exp or invexp

x = [0:0.05:1];

_exp = 0.26 .* exp (2.3.*x) + 1.05;
_rev_exp = 0.26 .* exp (2.3.*(1-x)) + 1.05;
_deriv_exp = 0.26 .* 2.3 .* exp (2.3.*x);
## Scaling for _exp to give inverse that starts and ends at the same point
A = 4.7727;
B =  -0.00014;
_inv_exp = A .* (1./_exp) + B;

## May also use a linear function
_lin = 1.31 + 2.333 .* x;

if strcmp (ligand_expression, 'exp')
  threshold = 1.1
else
  threshold = 1.1 * (0.26 .* exp(2.3) + 1.05)
end

## Debugging the curves
figure(7); clf;
hold on;
plot (x, _exp);
plot (x, flip(1./_exp));
plot (x, flip(_inv_exp));
plot (x, _lin);
legend('exp', '1/exp (flipped LR)','scaled 1/exp (flipped LR)', 'linear')
figure(1);

%% The J effect
_r0 = flip(_exp);
## Either this, making no strong assumptions about ligand expression:
if strcmp (ligand_expression, 'exp')
  _l0 = _exp;
elseif strcmp (ligand_expression, 'lin')
  _l0 = _lin;
else
  _l0 = flip(_inv_exp);
end

subplot(2,4,1);
plot (x, _r0, 'linestyle', '--'); % exp for receptors on retina
hold on;
plot (x, _l0, 'linestyle', '--'); % ligands
%%plot (x, _r0 .* _l0, 'color', colour.orchid2, 'linestyle', ':');
plot (x, _r0(1) .* _l0, 'r-');
plot (x(3:7), _r0(5) .* _l0(3:7), 'b');
plot (x(7:11), _r0(9) .* _l0(7:11), 'g');
plot (x(11:15), _r0(13) .* _l0(11:15), 'm');
plot (x(15:19), _r0(17) .* _l0(15:19), 'c');

plot ([x(1),x(1)], [_r0(1),_r0(1).*_l0(1)], 'ro-');
plot ([x(5),x(5)], [_r0(5),_r0(5).*_l0(5)], 'bo-');
plot ([x(9),x(9)], [_r0(9),_r0(9).*_l0(9)], 'go-');
plot ([x(13),x(13)], [_r0(13),_r0(13).*_l0(13)], 'mo-');
plot ([x(17),x(17)], [_r0(17),_r0(17).*_l0(17)], 'co-');

plot (x, _r0(21) .* _l0, 'k');
plot ([x(21),x(21)], [_r0(21),_r0(21).*_l0(21)], 'ko-');
lg1 = legend (['r0';'l0';'r0[0] x l0 interaction (competitive)';'r0[.2] x l0 interaction (competitive)'],'Location','North');
strlbl = ['Temporal ------------------------> Nasal'];
xlabel(strlbl)
ylabel('Expression/Interaction')
title('Mass-action r0/l0 expression/interaction')

_r2 = _exp; % opposing ligands
if strcmp (ligand_expression, 'exp')
  _l2 = flip(_exp); % opposing receptors
elseif strcmp (ligand_expression, 'lin')
  _l2 = flip(_lin);
else
  _l2 = _inv_exp;
end
subplot(2,4,2);
plot (x, _r2, 'linestyle', ':');
hold on;
plot (x, _l2, 'linestyle', ':'); % ligands
%%plot (x, _r2 .* _l2, 'color', colour.orchid2, 'linestyle', ':');
%%plot (x(1:3), _r2(1) .* _l2(1:3), 'r-');
plot (x, _r2(1) .* _l2, 'r');
plot ([x(1),x(1)], [_r2(2),_r2(1).*_l2(1)], 'ro-');
plot (x(3:7), _r2(5) .* _l2(3:7), 'b');
plot ([x(5),x(5)], [_r2(5),_r2(5).*_l2(5)], 'bo-');
plot (x(7:11), _r2(9) .* _l2(7:11), 'g');
plot ([x(9),x(9)], [_r2(9),_r2(9).*_l2(9)], 'go-');
plot (x(11:15), _r2(13) .* _l2(11:15), 'm');
plot ([x(13),x(13)], [_r2(13),_r2(13).*_l2(13)], 'mo-');
plot (x(15:19), _r2(17) .* _l2(15:19), 'c');
%%plot (x, _r2(17) .* _l2, 'k');
plot ([x(17),x(17)], [_r2(17),_r2(17).*_l2(17)], 'co-');

plot (x, _r2(21) .* _l2, 'k');
plot ([x(21),x(21)], [_r2(21),_r2(21).*_l2(21)], 'ko-');

lg2 = legend (['rcpt2';'lgnd2'],'Location','North');
strlbl = ['Temporal ------------------------> Nasal'];
##ylim([0,20])
xlabel(strlbl)
ylabel('Expression/Interaction')
title('Mass-action r2/l2 expression/interaction')

%% Now plot combined repulsion due to both receptor/ligand pairs

subplot(2,4,3);
hold on;
plot (x, _r2(1) .* _l2 + _r0(1) .* _l0, 'r-');
plot (x, _r2(5) .* _l2 + _r0(5) .* _l0, 'b-');
plot (x, _r2(9) .* _l2 + _r0(9) .* _l0, 'g-');
plot (x, _r2(13) .* _l2 + _r0(13) .* _l0, 'm-');
plot (x, _r2(17) .* _l2 + _r0(17) .* _l0, 'c-');
plot (x, _r2(21) .* _l2 + _r0(21) .* _l0, 'k-');
plot ([x(1),x(1)], [0,_r2(1).*_l2(1)+_r0(1).*_l0(1)], 'ro-');
plot ([x(5),x(5)], [0,_r2(5).*_l2(5)+_r0(5).*_l0(5)], 'bo-');
plot ([x(9),x(9)], [0,_r2(9).*_l2(9)+_r0(9).*_l0(9)], 'go-');
plot ([x(13),x(13)], [0,_r2(13).*_l2(13)+_r0(13).*_l0(13)], 'mo-');
plot ([x(17),x(17)], [0,_r2(17).*_l2(17)+_r0(17).*_l0(17)], 'co-');
plot ([x(21),x(21)], [0,_r2(21).*_l2(21)+_r0(21).*_l0(21)], 'ko-');
lg3 = legend (['r0[0]xl0 + r2[0]xl2 combined interaction (competitive)';'r0[.2]xl0 + r2[.2]xl2 combined interaction (competitive)';'etc'],'Location','North');
strlbl = ['Temporal ------------------------> Nasal'];
##ylim([0,20])
xlabel(strlbl)
ylabel('Expression/Interaction')
title('Summed mass-action interaction')

subplot(2,4,4);
hold on;
###
plot (x, _r0, 'linestyle', '--'); % exp for receptors on retina
hold on;
plot (x, _l0, 'linestyle', '--'); % ligands
plot (x, _r0(1) .* _l0, 'r-');
plot (x, _r0(5) .* _l0, 'b');
plot (x(7:11), _r0(9) .* _l0(7:11), 'g');
plot (x, _r0(13) .* _l0, 'm');
plot (x(15:19), _r0(17) .* _l0(15:19), 'c');

plot ([x(1),x(1)], [_r0(1),_r0(1).*_l0(1)], 'ro-');
plot ([x(5),x(5)], [_r0(5),_r0(5).*_l0(5)], 'bo-');
plot ([x(9),x(9)], [_r0(9),_r0(9).*_l0(9)], 'go-');
plot ([x(13),x(13)], [_r0(13),_r0(13).*_l0(13)], 'mo-');
plot ([x(17),x(17)], [_r0(17),_r0(17).*_l0(17)], 'co-');

plot (x, _r0(21) .* _l0, 'k');
plot ([x(21),x(21)], [_r0(21),_r0(21).*_l0(21)], 'ko-');
###
plot (x, _r2, 'linestyle', ':');
hold on;
plot (x, _l2, 'linestyle', ':'); % ligands
%%plot (x, _r2 .* _l2, 'color', colour.orchid2, 'linestyle', ':');
%%plot (x(1:3), _r2(1) .* _l2(1:3), 'r-');
plot (x, _r2(1) .* _l2, 'r');
plot ([x(1),x(1)], [_r2(2),_r2(1).*_l2(1)], 'ro-');
plot (x, _r2(5) .* _l2, 'b');
plot ([x(5),x(5)], [_r2(5),_r2(5).*_l2(5)], 'bo-');
plot (x(7:11), _r2(9) .* _l2(7:11), 'g');
plot ([x(9),x(9)], [_r2(9),_r2(9).*_l2(9)], 'go-');
plot (x, _r2(13) .* _l2, 'm');
plot ([x(13),x(13)], [_r2(13),_r2(13).*_l2(13)], 'mo-');
plot (x(15:19), _r2(17) .* _l2(15:19), 'c');
%%plot (x, _r2(17) .* _l2, 'k');
plot ([x(17),x(17)], [_r2(17),_r2(17).*_l2(17)], 'co-');

plot (x, _r2(21) .* _l2, 'k');
plot ([x(21),x(21)], [_r2(21),_r2(21).*_l2(21)], 'ko-');

plot ([0,1],[threshold,threshold],'k:')

###
strlbl = ['Temporal ------------------------> Nasal'];
xlabel(strlbl)
ylabel('Expression/Interaction')
title('Combined mass-action interaction')


## Now the retina->tectum interactions
subplot(2,4,5);
hold on;
_r0 = flip(_exp);
_L0 = _deriv_exp;
plot (x, _r0, 'linestyle', '--'); % exp for receptors on retina
plot (x, _L0, 'linestyle', '--'); % ligands
%%plot (x, _r0 .* _L0, 'color', colour.orchid2, 'linestyle', ':');
plot (x, _r0(1) .* _L0, 'r-');
plot (x(1:7), _r0(5) .* _L0(1:7), 'b');
plot (x(1:11), _r0(9) .* _L0(1:11), 'g');
plot (x(1:15), _r0(13) .* _L0(1:15), 'm');
plot (x(1:19), _r0(17) .* _L0(1:19), 'c');

plot ([x(1),x(1)], [_r0(1),_r0(1).*_L0(1)], 'ro-');
plot ([x(5),x(5)], [_r0(5),_r0(5).*_L0(5)], 'bo-');
plot ([x(9),x(9)], [_r0(9),_r0(9).*_L0(9)], 'go-');
plot ([x(13),x(13)], [_r0(13),_r0(13).*_L0(13)], 'mo-');
plot ([x(17),x(17)], [_r0(17),_r0(17).*_L0(17)], 'co-');

plot (x, _r0(21) .* _L0, 'k');
plot ([x(21),x(21)], [_r0(21),_r0(21).*_L0(21)], 'ko-');
lg4 = legend (['r0';'gradL0';'r0[0] x gradL0 interaction (tectal)';'r0[.2] x gradL0 interaction (tectal)'],'Location','North');
strlbl = ['Temporal/Rostral ------------------------> Nasal/Caudal'];
xlabel(strlbl)
ylabel('Gradient interaction')
title ('Receptor-tectal ligand interaction r0/L0')

subplot(2,4,6);
hold on;
_r2 = _exp;
_L2 = -flip(_deriv_exp);
plot (x, _r2, 'linestyle', '--'); % exp for receptors on retina
plot (x, _L2, 'linestyle', '--'); % ligands

plot (x, _r2(1) .* _L2, 'r-');
plot (x(1:7), _r2(5) .* _L2(1:7), 'b');
plot (x(1:11), _r2(9) .* _L2(1:11), 'g');
plot (x(1:15), _r2(13) .* _L2(1:15), 'm');
plot (x(1:19), _r2(17) .* _L2(1:19), 'c');

plot ([x(1),x(1)], [_r2(1),_r2(1).*_L2(1)], 'ro-');
plot ([x(5),x(5)], [_r2(5),_r2(5).*_L2(5)], 'bo-');
plot ([x(9),x(9)], [_r2(9),_r2(9).*_L2(9)], 'go-');
plot ([x(13),x(13)], [_r2(13),_r2(13).*_L2(13)], 'mo-');
plot ([x(17),x(17)], [_r2(17),_r2(17).*_L2(17)], 'co-');

plot (x, _r2(21) .* _L2, 'k');
plot ([x(21),x(21)], [_r2(21),_r2(21).*_L2(21)], 'ko-');
lg5 = legend (['r2';'gradL2';'r2[0] x grad L2 interaction (tectal)';'r2[.2] x grad L2 interaction (tectal)'],'Location','Southeast');
strlbl = ['Temporal/Rostral ------------------------> Nasal/Caudal'];
xlabel(strlbl)
ylabel('Gradient interaction')
title ('Receptor-tectal ligand interaction r2/L2')


subplot(2,4,7);
hold on;
plot (x, _r2(1) .* _L2 + _r0(1) .* _L0, 'r-');
plot (x, _r2(5) .* _L2 + _r0(5) .* _L0, 'b-');
plot (x, _r2(9) .* _L2 + _r0(9) .* _L0, 'g-');
plot (x, _r2(13) .* _L2 + _r0(13) .* _L0, 'm-');
plot (x, _r2(17) .* _L2 + _r0(17) .* _L0, 'c-');
plot (x, _r2(21) .* _L2 + _r0(21) .* _L0, 'k-');
plot ([x(1),x(1)], [0,_r2(1).*_L2(1)+_r0(1).*_L0(1)], 'ro-');
plot ([x(5),x(5)], [0,_r2(5).*_L2(5)+_r0(5).*_L0(5)], 'bo-');
plot ([x(9),x(9)], [0,_r2(9).*_L2(9)+_r0(9).*_L0(9)], 'go-');
plot ([x(13),x(13)], [0,_r2(13).*_L2(13)+_r0(13).*_L0(13)], 'mo-');
plot ([x(17),x(17)], [0,_r2(17).*_L2(17)+_r0(17).*_L0(17)], 'co-');
plot ([x(21),x(21)], [0,_r2(21).*_L2(21)+_r0(21).*_L0(21)], 'ko-');
plot ([0,1],[0,0],'k:')
lg6 = legend (['r0[0]xL0 + r2[0]xL2 combined interaction (tectal)';'r0[.2]xL0 + r2[.2]xL2 combined interaction (tectal)';'etc'],'Location','North');
strlbl = ['Temporal/Rostral ------------------------> Nasal/Caudal'];
xlabel(strlbl)
ylabel('Gradient interaction')


subplot(2,4,8);
hold on;
xx = [0:0.05:0.5];
_L2trunc = _L2(1:11);
_L0trunc = _L0(1:11);

plot (xx, _r2(1) .* _L2trunc + _r0(1) .* _L0trunc, 'r-');
plot (xx, _r2(5) .* _L2trunc + _r0(5) .* _L0trunc, 'b-');
plot (xx, _r2(9) .* _L2trunc + _r0(9) .* _L0trunc, 'g-');
plot (xx, _r2(13) .* _L2trunc + _r0(13) .* _L0trunc, 'm-');
plot (xx, _r2(17) .* _L2trunc + _r0(17) .* _L0trunc, 'c-');
plot (xx, _r2(21) .* _L2trunc + _r0(21) .* _L0trunc, 'k-');
plot ([xx(1),xx(1)], [0,_r2(1).*_L2(1)+_r0(1).*_L0(1)], 'ro-');
plot ([xx(3),xx(3)], [0,_r2(5).*_L2(2)+_r0(5).*_L0(3)], 'bo-');
plot ([xx(5),xx(5)], [0,_r2(9).*_L2(5)+_r0(9).*_L0(5)], 'go-');
plot ([xx(7),xx(7)], [0,_r2(13).*_L2(7)+_r0(13).*_L0(7)], 'mo-');
plot ([xx(9),xx(9)], [0,_r2(17).*_L2(9)+_r0(17).*_L0(9)], 'co-');
plot ([xx(11),xx(11)], [0,_r2(21).*_L2(11)+_r0(21).*_L0(11)], 'ko-');
plot ([0,0.5],[0,0],'k:')
lg7 = legend (['r0[0]xL0 + r2[0]xL2 combined interaction (tectal)';'r0[.2]xL0 + r2[.2]xL2 combined interaction (tectal)';'etc'],'Location','South');
strlbl = ['Temporal/Rostral -----------------------------> Nasal/Tectal Midline'];
xlabel(strlbl)





#figure(100)
#clf;
#plot (x, _exp)
#hold on
#plot (x, _deriv_exp, ':')
