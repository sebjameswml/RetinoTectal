% Reber et al 2004, Box 1 Equations

x = [0:0.1:100];

% Equation 4 wildtype: +/+
plus_plus = 0.26 .* exp (0.023.*x) + 1.05;

% Equation 5. heterozygous EphA3 knockin: ki/+
%ki_plus =  0.26 .* exp (0.023.*x) + 1.98;
ki_plus =  plus_plus + (1.98-1.05);

% Equation 6. homozygous EphA3 knockin: ki/ki
ki_ki =   plus_plus + (2.91-1.05);

% Equation 7 - Local relative signalling ratio for ki/+ EphA3 knockin
%
% This is it - Reber et al say that the signal is ratio of knockin to
% wildtype expression. Really, this just says that at some point, the
% system can't discriminate between the levels. And is this true for
% retinal levels or is it due to the ligand levels on the tectum?
R_ki_plus = ki_plus ./ plus_plus;

% Eq. 8 - ratio for ki/ki EphA3 knockin
R_ki_ki = ki_ki ./ plus_plus;

figure(1); clf;
plot (x, plus_plus);
hold on;
plot (x, ki_plus);
plot (x, R_ki_plus);
legend (['+/+ (Eq 4)';'ki/+ (Eq 5)';'Signal ki/+ (Eq 7 = Eq 5/Eq 4)']);
