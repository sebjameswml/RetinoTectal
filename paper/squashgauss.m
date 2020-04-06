x = -3:0.01:3;

% Params
%m = 15       % the sharpness/gauss height combined multiplier
mu = 1       % mean of Gaussian - supplied by gamma_i,1 or gamma_i,2
sigma = 0.2  % Width of Gaussian; a parameter

figure (1); clf; hold on;

g = exp (-0.5 .* power (((x-mu)./sigma), 2));
f = 2 .* (-0.5+(1./(1 + exp(-g.*m))));

plot (x,g,'-b');
plot (x,1-f,'-r');
lg = legend('Gaussian','Inverted squashed Gaussian')
set (lg, 'Location', 'NorthWest')