x = -3:0.01:3;

% Params
m = 5        % multiplier on Gaussian; make sure it's >1; a parameter
mu = 1       % mean of Gaussian - supplied by gamma_i,1 or gamma_i,2
sigma = 0.2  % Width of Gaussian; a parameter
sh = 3       % sharpness of sigmoid; a parameter

% (m * sh) is really one parameter because if you are only
% interested in f(x), and you write it out fully, then m simply
% multiplies s. Thus, there are 2 tunable parameters and the
% positional parameter, mu.

figure (1); clf; hold on;

g = m .* exp (-0.5 .* power (((x-mu)./sigma), 2));
f = 2 .* (-0.5+(1./(1 + exp(-g.*sh))));

plot (x,g,'-b');
plot (x,1-f,'-r');
lg = legend('Gaussian','Inverted squashed Gaussian')
set (lg, 'Location', 'NorthWest')