
% First manually load data from the run:

% E.g. load log/anneal1/anneal1_ee_GI_20210928_160433.h5

pa = param_hist_accepted;
pr = param_hist_rejected;
pb = [pa, pr];

fpa = f_param_hist_accepted;
fpr = f_param_hist_rejected;
fpb = [fpa, fpr];

lim_maxes = max(pb');
lim_mins = min(pb');

fpanorm = fpa./max(fpa);
fprnorm = fpr./max(fpr);
fpbnorm = fpb./max(fpb);

fnumstart = 0;

% Scatter plot with larger size correlating with larger objective
% value and 'hotter' colour correlating with smaller (i.e. better)
% objective value. So small hot points are better.
figure(fnumstart+1);
scatter3(pa(1,:), pa(2,:), pa(3,:), 2000.*fpanorm, (1-fpanorm),
         "filled");
xlim([lim_mins(1), lim_maxes(1)]);
ylim([lim_mins(2), lim_maxes(2)]);
zlim([lim_mins(3), lim_maxes(3)]);
title ('accepted');
xlabel(param_name_1')
ylabel(param_name_2')
zlabel(param_name_3')

figure(fnumstart+2);
scatter3(pr(1,:), pr(2,:), pr(3,:), 2000.*fprnorm, (1-fprnorm), "filled");
xlim([lim_mins(1), lim_maxes(1)]);
ylim([lim_mins(2), lim_maxes(2)]);
zlim([lim_mins(3), lim_maxes(3)]);
title ('rejected');
xlabel(param_name_1')
ylabel(param_name_2')
zlabel(param_name_3')

figure(fnumstart+3);
scatter3(pb(1,:), pb(2,:), pb(3,:), 2000.*fpbnorm, (1-fpbnorm), "filled");
xlim([lim_mins(1), lim_maxes(1)]);
ylim([lim_mins(2), lim_maxes(2)]);
zlim([lim_mins(3), lim_maxes(3)]);
title ('both');
xlabel(param_name_1')
ylabel(param_name_2')
zlabel(param_name_3')
