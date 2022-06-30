##
## Epha3/Epha4
##
function examineRetEph (knockin_ephax = 0, knockdown_epha4 = 0, fn)

  h_f = figure(fn);
  h_f_pos = get(h_f, 'Position');
  ##                                            w     h
  set(h_f, 'Position', [h_f_pos(1), h_f_pos(2), 2000, 600]);
  clf;

  x = [0:0.05:1];

  w_EphCommon = 0.25;

  _exp = 0.26 .* exp (2.3.*x) + 1.05; % T exponential_expression (const T& x)

  ## Start with a certain amount of ephrinA
  ephrinA = fliplr(_exp);


  ## Even expression of EphA4
  _epha4 = ones(1,length(ephrinA)) .* 2.0 - knockdown_epha4;

  ## Binding affinity for EphA4
  w_EphA4 = 0.25;

  ## Let EphA4 get first dibs on ephrinA
  _p_epha4 = (w_EphA4 .* _epha4 .* ephrinA);

  ## Remaining epha4 could interact
  EphA4_free = _epha4 - _p_epha4;

  ephrinA_free = ephrinA - _p_epha4;


  ## And some expression of EphA3/x whatever
  EphAx = _exp + knockin_ephax;

  ## Binding affinity for EphAx
  w_EphAx = w_EphCommon;

  ## Some ephrinA binds to EphAx. Really, need a first order model that takes into account amount of EphAx and amount of ephrinA
  ephrinA_binds_EphAx = w_EphAx .* EphAx .* ephrinA_free;

  ## Remaining free receptors can interact with tectal ephrins
  EphAx_free = EphAx - ephrinA_binds_EphAx;

  ephrinA_free2 = ephrinA_free - ephrinA_binds_EphAx;


  subplot(1,5,1);
  hold on;
  plot (x, EphAx, 'linestyle', ':', 'color', 'blue'); % exp for receptors on retina
  plot (x, EphAx_free, 'linestyle', '-', 'color', 'blue');
  plot (x, ephrinA_binds_EphAx, 'linestyle', '--', 'color', 'blue');
  lg1 = legend (['EphAx';'EphAx free';'EphAx bound'], 'Location','North');
  strlbl = ['Nasal ------------------------> Temporal'];
  ylim([0,5])
  str = sprintf('EphAx knockin %.1f', knockin_ephax);
  text(0.6, 0.9, str)
  str = sprintf('EphA4 knockdown %.1f', knockdown_epha4);
  text(0.6, 0.5, str)
  xlabel(strlbl)
  ylabel('Expression/Interaction')

  subplot(1,5,2);
  hold on;
  plot (x, ephrinA, 'linestyle', ':', 'color', 'black'); % ligands
  plot (x, ephrinA_free, 'linestyle', '-', 'color', 'black'); % ligands
  plot (x, ephrinA_binds_EphAx, 'linestyle', '--', 'color', 'blue'); % ligands
  plot (x, _p_epha4, 'linestyle', '--', 'color', 'red'); % ligands
  lg2 = legend (['ephrinA';'ephrinA free for EphAx';'ephrinA bound to EphA3';'ephrinA bound to EphA4'], 'Location','North');
  ylim([0,5])
  xlabel(strlbl)
  ylabel('Expression/Interaction')

  subplot(1,5,3);
  hold on;
  plot (x, _epha4, 'linestyle', ':', 'color', 'red');
  plot (x, EphA4_free, 'linestyle', '-', 'color', 'red'); #not phosphorylized = available
  plot (x, _p_epha4, 'linestyle', '--', 'color', 'red');
  lg1 = legend (['EphA4';'EphA4 free';'EphA4 bound'], 'Location','North');
  ylim([0,5])
  xlabel(strlbl)
  ylabel('Expression/Interaction')

  axpow = 1;
  ###############################
  subplot(1,5,4);
  hold on;

  # Re-write EphA4_free as quadratic?
  #EphA4_free = 1.31 + 2.3333 .* x .* x; # This is what I have in the sim right now
  #EphA4_free = _p_epha4 # That's really _bound_ EphA4

  ax = EphAx_free./EphA4_free;

  plot (x, 0.01 .* power(ax, axpow))
  plot (x, EphAx_free .* 0.01 .* power(ax, axpow))

  lg1 = legend (['(Ax/A4)^k';'Ax * (Ax/A4)^k'], 'Location','North');
  xlabel(strlbl)
  ylabel('Expression/Interaction')

  ################
  subplot(1,5,5);
  hold on;

  plot (x, EphAx_free)
  plot (x, EphA4_free)
  plot (x, EphAx_free + EphAx_free .* 0.1 .* power(ax, axpow))
  plot (x, ax)

  lg1 = legend (['Ax';'A4';'Ax(1 + (Ax/A4)^k)';'Ax/A4'], 'Location','North');
  xlabel(strlbl)
  ylabel('Expression/Interaction')

end