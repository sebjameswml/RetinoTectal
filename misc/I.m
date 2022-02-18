% Equations used in my modelling. See tissue.h and
% guidingtissue::exponential_expression and similar

% This models the S&G-style receptor-receptor interaction based on receptor ratios

## Use sebcolour.m (must be in path)
sebcolour

x = [0:0.05:1];

_exp = 0.26 .* exp (2.3.*x) + 1.05;
_exp2 = 0.26 .* exp (1.1.*x) + 1.05;
_log = 2.32 + 1.29 .* log (2.3 .* (x+0.2));
_quad = 1.31 + 2.333 .* x .* x;
_lin = 1.31 + 2.333 .* x;


%% The I effect
_r0 = flip(_exp); % One receptor type.
_r2 = _exp;       % A second receptor type.
figure(10); clf;
plot (x, _r0, 'linestyle', '--'); % exp for receptors on retina
hold on;
plot (x, _r2, 'linestyle', '--', 'color', colour.gray90); % ligands

plot ([0,1], [1.1,1.1], 'k:');

plot (x, _r0(1) ./ _r0, 'r-');
plot (x, _r0(5) ./ _r0, 'b-');
plot (x, _r0(9) ./ _r0, 'g-');
plot (x, _r0(13) ./ _r0, 'm-');
plot (x, _r0(17) ./ _r0, 'c-');
plot (x, _r0(21) ./ _r0, 'k-');


## Vertical lines
plot ([x(1),x(1)], [_r0(1), _r0(1)./_r0(1)], 'ro-');
plot ([x(5),x(5)], [_r0(5), _r0(5)./_r0(5)], 'bo-');
plot ([x(9),x(9)], [_r0(9), _r0(9)./_r0(9)], 'go-');
plot ([x(13),x(13)], [_r0(13), _r0(13)./_r0(13)], 'mo-');
plot ([x(17),x(17)], [_r0(17), _r0(17)./_r0(17)], 'co-');
plot ([x(21),x(21)], [_r0(21), _r0(21)./_r0(21)], 'ko-');

lg = legend (['r0';'r2';'1.1 threshold';
              'r0[b=0] / r0[k] (receptor ratio)';
              'r0[b=.2] / r0[k] (receptor ratio)';
              'r0[b=.4] / r0[k] (receptor ratio)'
             ],
             'Location','North');



strlbl = ['Temporal -------------------------------------> Nasal'];
figure(10)
##ylim([0,20])
xlabel(strlbl)


figure(11); clf;
plot (x, _r0, 'linestyle', '--', 'color', colour.gray90); % exp for receptors on retina
hold on;
plot (x, _r2, 'linestyle', '--'); % ligands

plot ([0,1], [1.1,1.1], 'k:');

plot (x, _r2(1) ./ _r2, 'r-');
plot (x, _r2(5) ./ _r2, 'b-');
plot (x, _r2(9) ./ _r2, 'g-');
plot (x, _r2(13) ./ _r2, 'm-');
plot (x, _r2(17) ./ _r2, 'c-');
plot (x, _r2(21) ./ _r2, 'k-');


## Vertical lines
plot ([x(1),x(1)], [_r2(1), _r2(1)./_r2(1)], 'ro-');
plot ([x(5),x(5)], [_r2(5), _r2(5)./_r2(5)], 'bo-');
plot ([x(9),x(9)], [_r2(9), _r2(9)./_r2(9)], 'go-');
plot ([x(13),x(13)], [_r2(13), _r2(13)./_r2(13)], 'mo-');
plot ([x(17),x(17)], [_r2(17), _r2(17)./_r2(17)], 'co-');
plot ([x(21),x(21)], [_r2(21), _r2(21)./_r2(21)], 'ko-');

lg = legend (['r0';'r2';'1.1 threshold';
              'r0[b=0] / r0[k] (receptor ratio)';
              'r0[b=.2] / r0[k] (receptor ratio)';
              'r0[b=.4] / r0[k] (receptor ratio)'
             ],
             'Location','North');

figure(12); clf;

plot (x, _r0, 'linestyle', '--'); % exp for receptors on retina
hold on;
plot (x, _r2, 'linestyle', '--'); % ligands

plot ([0,1], [1.1,1.1], 'k:');

plot (x, _r2(1) ./ _r2, 'r-');
plot (x, _r2(5) ./ _r2, 'b-');
plot (x, _r2(9) ./ _r2, 'g-');
plot (x, _r2(13) ./ _r2, 'm-');
plot (x, _r2(17) ./ _r2, 'c-');
plot (x, _r2(21) ./ _r2, 'k-');

## Vertical lines
plot ([x(1),x(1)], [_r2(1), _r2(1)./_r2(1)], 'ro-');
plot ([x(5),x(5)], [_r2(5), _r2(5)./_r2(5)], 'bo-');
plot ([x(9),x(9)], [_r2(9), _r2(9)./_r2(9)], 'go-');
plot ([x(13),x(13)], [_r2(13), _r2(13)./_r2(13)], 'mo-');
plot ([x(17),x(17)], [_r2(17), _r2(17)./_r2(17)], 'co-');
plot ([x(21),x(21)], [_r2(21), _r2(21)./_r2(21)], 'ko-');

plot (x, _r0(1) ./ _r0, 'r-');
plot (x, _r0(5) ./ _r0, 'b-');
plot (x, _r0(9) ./ _r0, 'g-');
plot (x, _r0(13) ./ _r0, 'm-');
plot (x, _r0(17) ./ _r0, 'c-');
plot (x, _r0(21) ./ _r0, 'k-');


## Vertical lines
plot ([x(1),x(1)], [_r0(1), _r0(1)./_r0(1)], 'ro-');
plot ([x(5),x(5)], [_r0(5), _r0(5)./_r0(5)], 'bo-');
plot ([x(9),x(9)], [_r0(9), _r0(9)./_r0(9)], 'go-');
plot ([x(13),x(13)], [_r0(13), _r0(13)./_r0(13)], 'mo-');
plot ([x(17),x(17)], [_r0(17), _r0(17)./_r0(17)], 'co-');
plot ([x(21),x(21)], [_r0(21), _r0(21)./_r0(21)], 'ko-');

title ('Any pair of branches with r0 and r2 expression will exhibit competitition!')
