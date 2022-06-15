##
## Epha3/Epha4
##
figure(100); clf;

x = [0:0.05:1];
_exp = 0.26 .* exp (2.3.*x) + 1.05; % T exponential_expression (const T& x)
ephrinA = fliplr(_exp);
_r0 = _exp;

_epha4 = ones(1,length(ephrinA)) .* 2.0;
_p_epha4 = _epha4 - (0.5 .* ephrinA);

subplot(2,2,1);
plot (x, _r0, 'linestyle', '-', 'color', 'blue'); % exp for receptors on retina
hold on;
plot (x, ephrinA, 'linestyle', ':', 'color', 'black'); % ligands
plot (x, _epha4, 'linestyle', '--', 'color', 'red');
plot (x, _p_epha4, 'linestyle', ':', 'color', 'red');
plot (x, _epha4 - _p_epha4, 'linestyle', '-', 'color', 'red');
lg1 = legend (['r0 (EphAx)';'ephrinA (ret)';'EphA4';'Phos. EphA4';'Avail. EphA4'], 'Location','North');
strlbl = ['Nasal ------------------------> Temporal'];
xlabel(strlbl)
ylabel('Expression/Interaction')

# Ki/Kd
subplot(2,2,2);
_r0_ki = _r0 + 1.0;
_epha4_kd = _epha4 * 0.74;
_p_epha4 = _epha4_kd - (0.5 .* ephrinA);

plot (x, _r0_ki, 'linestyle', '-', 'color', 'blue'); % exp for receptors on retina
hold on;
plot (x, ephrinA, 'linestyle', ':', 'color', 'black'); % ligands
plot (x, _epha4_kd, 'linestyle', '--', 'color', 'red');
plot (x, _p_epha4, 'linestyle', ':', 'color', 'red');
plot (x, _epha4_kd - _p_epha4, 'linestyle', '-', 'color', 'red');
lg1 = legend (['r0 (EphAx)';'ephrinA (ret)';'EphA4';'Phos. EphA4';'Avail. EphA4'], 'Location','North');
strlbl = ['Nasal ------------------------> Temporal'];
xlabel(strlbl)
ylabel('Expression/Interaction')
title('Knock in/out');
