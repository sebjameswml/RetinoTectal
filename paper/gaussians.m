x = -3:0.01:3;

% Params
m = 5        % multiplier on Gaussian; make sure it's >1
sigma = 0.2   % Width of Gaussian
sh = 3      % sharpness of sigmoid

figure (1); clf; hold on;

g = m .* exp (-0.5 .* power (((x-mu)./sigma), 2));
plot (x,g);
gg = 2 .* (-0.5+(1./(1 + exp(-g.*sh))));
plot (x,1-gg)
