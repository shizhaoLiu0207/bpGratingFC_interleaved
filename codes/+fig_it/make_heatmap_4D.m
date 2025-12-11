function make_heatmap_4D(feature_target, cardinal_delta,cardinal_prior, oblique_delta,oblique_prior, plotOptions)


subplot(2,3,1)
p1 = cardinal_delta;
p2 = cardinal_prior;
p1_name = '\color{red}\delta_{cardinal}';
p2_name = '\color{red}prior_{cardinal}';

fig_it.make_heatmap(p1,p2,feature_target,p1_name,p2_name, plotOptions)

set(gca,'fontsize',plotOptions.ftsize)


subplot(2,3,2)
p1 = cardinal_delta;
p2 = oblique_delta;
p1_name = '\color{red}\delta_{cardinal}';
p2_name = '\color{blue}\delta_{oblique}';
fig_it.make_heatmap(p1,p2,feature_target,p1_name,p2_name, plotOptions)
set(gca,'fontsize',plotOptions.ftsize)

subplot(2,3,3)
p1 = cardinal_delta;
p2 = oblique_prior;
p1_name = '\color{red}\delta_{cardinal}';
p2_name = '\color{blue}prior_{oblique}';
fig_it.make_heatmap(p1,p2,feature_target,p1_name,p2_name, plotOptions)
set(gca,'fontsize',plotOptions.ftsize)

subplot(2,3,4)
p1 = cardinal_prior;
p2 = oblique_delta;
p1_name = '\color{red}prior_{cardinal}';
p2_name = '\color{blue}\delta_{oblique}';
fig_it.make_heatmap(p1,p2,feature_target,p1_name,p2_name, plotOptions)
set(gca,'fontsize',plotOptions.ftsize)

subplot(2,3,5)
p1 = cardinal_prior;
p2 = oblique_prior;
p1_name = '\color{red}prior_{cardinal}';
p2_name = '\color{blue}prior_{oblique}';
fig_it.make_heatmap(p1,p2,feature_target,p1_name,p2_name, plotOptions)
set(gca,'fontsize',plotOptions.ftsize)

subplot(2,3,6)
p1 = oblique_delta;
p2 = oblique_prior;
p1_name = '\color{blue}\delta_{oblique}';
p2_name = '\color{blue}prior_{oblique}';
fig_it.make_heatmap(p1,p2,feature_target,p1_name,p2_name, plotOptions)
set(gca,'fontsize',plotOptions.ftsize)
sgtitle(plotOptions.feature_title, 'fontweight','bold','fontsize',plotOptions.ftsize + 2,'interpreter','tex')
end