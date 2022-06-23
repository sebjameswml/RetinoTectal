## differential equations.
function xdot = de (x, t)

  xdot = zeros (3,1);

  ## Equation for EphAx
  gamma3 = 0.01;
  alpha3 = 0.0;
  xdot(1) = gamma3 * x(3) * x(1) - alpha3 * x(1);
  ## Equation for EphA4
  gamma4 = 0.02;
  alpha4 = 0.0;
  xdot(2) = gamma4 * x(3) * x(2) - alpha4 * x(2);
  ## Equation for ephrinA
  xdot(3) = -xdot(1)-xdot(2);

endfunction
