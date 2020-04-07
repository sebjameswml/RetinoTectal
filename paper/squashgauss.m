x = -3:0.01:3;

% Params
s = 15       % the sharpness/gauss height combined multiplier
mu = 1       % mean of Gaussian - supplied by gamma_i,1 or gamma_i,2
w = 0.2      % Width of Gaussian; a parameter

figure (1); clf; hold on;

g = exp (- power (((x-mu)./(2.*w)), 2));
f = 2 ./(1 + exp(-g.*s)) - 1;
finv = 1 - f;
plot (x,g,'-b');
plot (x,f,'-g');
plot (x,finv,'-r');
lg = legend('Gaussian','Squashed Gaussian','Inverted squashed Gaussian')
set (lg, 'Location', 'NorthWest')