% Algebraic example from Section 2.2 of Sensitivity Analysis for Periodic
% Orbits and Quasiperiodic Invariant Tori Using the Adjoint Method, by
% Dankowicz & Sieber

function inflection()
% Run with inactive continuation parameters 'av', 'ep', and 'ze':
prob = coco_prob();
prob = coco_add_func(prob, 'phi', @phi, @dphi, [], 'zero', ...
  'u0', [-0.49; 4.9; 4.9; 0; 5; 5; 1.01; 1; .1; .01]);
prob = coco_add_func(prob, 'psi', @psi, @dpsi, [], 'inactive', ...
  {'da', 'av', 'ep', 'ze'}, 'uidx', 1:10);
prob = coco_add_adjt(prob, 'phi');
prob = coco_add_adjt(prob, 'psi', {'e.da','e.av','e.ep','e.ze'}, ...
  'aidx', 1:10);
coco(prob, 'run', [], 1, {'da', 'e.da', 'e.av', 'e.ep', 'e.ze'}, ...
  {[], [0 1]})


% Run with inactive continuation parameters 'ep', 'ze', and 'e.da' set
% equal to 1:
prob = coco_prob();
chart = coco_read_solution('phi', 'run', 5, 'chart');
prob = coco_add_func(prob, 'phi', @phi, @dphi, [], 'zero', ...
  'u0', chart.x);
prob = coco_add_func(prob, 'psi', @psi, @dpsi, [], 'inactive', ...
  {'da', 'av', 'ep', 'ze'}, 'uidx', 1:10);
chart = coco_read_adjoint('phi', 'run', 5, 'chart');
prob = coco_add_adjt(prob, 'phi', 'l0', chart.x);
chart = coco_read_adjoint('psi', 'run', 5, 'chart');
prob = coco_add_adjt(prob, 'psi', {'e.da','e.av','e.ep','e.ze'}, ...
  'aidx', 1:10, 'l0', chart.x);
prob = coco_set(prob, 'cont', 'ItMX', 500, 'NPR', 100);
coco(prob, 'run', [], 1, {'da', 'av', 'e.av', 'e.ep', 'e.ze'}, ...
  {[], [0.5 2.5]})

end

%% zero and monitor functions

function [data, f] = phi(prob, data, u)

v = num2cell(u);
[a1, b1, c1, a2, b2, c2, o1, o2, ze, ep] = deal(v{:});
f = [c1^2-a1^2-b1^2; c2^2-a2^2-b2^2; o1-o2-ep;
  (1-o1^2)*a1+2*ze*o1*b1-1; (1-o1^2)*b1-2*ze*o1*a1;
  (1-o2^2)*a2+2*ze*o2*b2-1; (1-o2^2)*b2-2*ze*o2*a2];

end

function [data, J] = dphi(prob, data, u)

v = num2cell(u);
[a1, b1, c1, a2, b2, c2, o1, o2, ze, ep] = deal(v{:});
J = [-2*a1,-2*b1,2*c1,0,0,0,0,0,0,0;
  0,0,0,-2*a2,-2*b2,2*c2,0,0,0,0;
  0,0,0,0,0,0,1,-1,0,-1;
  1-o1^2,2*ze*o1,0,0,0,0,-2*o1*a1+2*ze*b1,0,2*o1*b1,0;
  -2*ze*o1,1-o1^2,0,0,0,0,-2*o1*b1-2*ze*a1,0,2*o1*a1,0;
  0,0,0,1-o2^2,2*ze*o2,0,0,-2*o2*a2+2*ze*b2,2*o2*b2,0;
  0,0,0,-2*ze*o2,1-o2^2,0,0,-2*o2*b2-2*ze*a2,2*o2*a2,0];

end

function [data, f] = psi(prob, data, u)

v = num2cell(u);
[a1, b1, c1, a2, b2, c2, o1, o2, ze, ep] = deal(v{:});
f = [c1-c2; (o1+o2)/2; ep; ze];

end

function [data, J] = dpsi(prob, data, u)

v = num2cell(u);
[a1, b1, c1, a2, b2, c2, o1, o2, ze, ep] = deal(v{:});
J = [0,0,1,0,0,-1,0,0,0,0;
  0,0,0,0,0,0,1/2,1/2,0,0;
  0,0,0,0,0,0,0,0,0,1;
  0,0,0,0,0,0,0,0,1,0];

end
