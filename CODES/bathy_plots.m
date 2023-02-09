%% Generate an example bathymetry plot to visualize all estimated bathymetries
% Pick two example hovers
r1 = 5
r2 = 29
textsz = 14;
Video_bathy=Video;
figure(1);clf
subplot(421)
rr = r1;
title(sprintf('Hover %i, MOP %i', Video_bathy(rr).flight, Video(rr).mop))
xmax = Video_bathy(rr).x10';xmin = Video_bathy(rr).x10'; zmax = (Video_bathy(rr).cbathy.z+Video_bathy(rr).cbathy.zerr)';zmin = (Video_bathy(rr).cbathy.z-Video_bathy(rr).cbathy.zerr)';
xmax(isnan(zmax)==1)=[]; zmax(isnan(zmax)==1)=[]; xmin(isnan(zmin)==1)=[]; zmin(isnan(zmin)==1)=[];
patch([xmax fliplr(xmin)], [zmax fliplr(zmin)],[0.6 0.6 0.6], 'EdgeColor', 'none');alpha 0.4
hold on
p(1)=plot(Video_bathy(rr).x10, Video_bathy(rr).survey.z, 'k', 'LineWidth', 4);
p(2)=plot(Video_bathy(rr).x10, Video_bathy(rr).cbathy.z, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
p(3)=plot(Video_bathy(rr).x10, Video_bathy(rr).cbathy.cbathy_hErr, 'r--', 'LineWidth', 3);
p(4)=plot(Video_bathy(rr).x10, Video_bathy(rr).cbathy.cbathy_gamma, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
p(5)=plot(Video_bathy(rr).x10, Video_bathy(rr).composite.cbathy_nlin,'--','Color', [0 0.4470 0.7410], 'LineWidth', 2);
lg = legend(p, 'Survey','Default', 'hErr < 0.5', 'Breaking', 'Nonlinear', 'Location', 'southwest', 'Box', 'off', 'FontSize', textsz-2);
xlabel('cBathy')
ylabel('Elevation (NAVD88 m)')
set(gca, 'FontSize', textsz)
grid on
xlim([0 400])
ylim([-8 2])
set(gca, 'box', 'on')

subplot(423)
p(1)=plot(Video_bathy(rr).x10, Video_bathy(rr).survey.z, 'k', 'LineWidth', 4);
hold on
p(2)=plot(Video_bathy(rr).x5, Video_bathy(rr).h_avg.lin,  'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
p(3)=plot(Video_bathy(rr).x5, Video_bathy(rr).h_avg.nlin, '--','Color', [0 0.4470 0.7410], 'LineWidth', 2);
p(4)=plot(Video_bathy(rr).x5, Video_bathy(rr).h_avg.bp, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 3);
lg = legend(p, 'Survey', 'Linear', 'Nonlinear', 'BP', 'Location', 'northeast', 'Box', 'off', 'FontSize', textsz-2, 'Location', 'southwest')
xlabel('Crest-tracking')
set(gca, 'FontSize', textsz)
ylabel('Elevation (NAVD88 m)')
grid on
xlim([0 400])
ylim([-8 2])

subplot(425)
p(1)=plot(Video_bathy(rr).x10, Video_bathy(rr).survey.z, 'k', 'LineWidth', 4);
hold on
p(2)=plot(Video_bathy(rr).x10, Video_bathy(rr).composite.cbathy_hErr, 'r--', 'LineWidth', 2);
p(3)=plot(Video_bathy(rr).x10, Video_bathy(rr).composite.cbathy_gamma,  'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
p(4)=plot(Video_bathy(rr).x10, Video_bathy(rr).composite.cbathy_nlin, '--','Color', [0 0.4470 0.7410], 'LineWidth', 2);
p(5)=plot(Video_bathy(rr).x10, Video_bathy(rr).composite.cbathyCT, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 3);
lg = legend(p, 'Survey', 'cB hErr < 0.5', 'cB Breaking', 'cB Nonlinear', 'cBathyCT', 'Location', 'southwest', 'Box', 'off', 'FontSize', textsz-2);
xlabel('Composite')
set(gca, 'FontSize', textsz)
ylabel('Elevation (NAVD88 m)')
grid on
xlim([0 400])
ylim([-8 2])

ax = subplot(427);
xticks([]); yticks([])
if Video_bathy(rr).flight < 10
    [fi,fcmap] = imread(fullfile(data_dir, 'timestacks','data', [char(string(Video_bathy(rr).date)) '_' char(Video_bathy(rr).location) '_0' char(string(Video_bathy(rr).flight)) '_' char(string(Video_bathy(rr).mop)) '.png']));
else
    [fi,fcmap] = imread(fullfile(data_dir, 'timestacks','data', [char(string(Video_bathy(rr).date)) '_' char(Video_bathy(rr).location) '_' char(string(Video_bathy(rr).flight)) '_' char(string(Video_bathy(rr).mop)) '.png']));
end 

ax2 = axes('Position',[ax.Position(1) ax.Position(2) ax.Position(3) ax.Position(4)/2]);
image(ax2,Video_bathy(rr).x10, 1,imrotate(cast(fi, 'uint8'),-90))
set(ax2,'ytick', [])
xlabel('Cross-shore Distance (m)', 'FontSize', textsz-4);
xlim([0 400])
   
ax3 = axes('Position',[ax.Position(1) ax.Position(2)+ax.Position(4)/2 ax.Position(3) ax.Position(4)/2]);
plot(ax3, Video_bathy(rr).x10, Video_bathy(rr).gamma_mean,'Color', [0.4660 0.6740 0.1880], 'LineWidth', 3)
hold on
plot(ax3, [100 100], [0 0.5], 'Color', [0.8 0.8 0.8])
plot(ax3, [200 200], [0 0.5], 'Color', [0.8 0.8 0.8])
plot(ax3, [300 300], [0 0.5], 'Color', [0.8 0.8 0.8])
plot(ax3, [400 400], [0 0.5], 'Color', [0.8 0.8 0.8])
ylim(ax3, [0 0.42])
set(ax3,'xtick', [])
xlim(ax3, [0 400])
set(ax2, 'FontSize', textsz)
set(ax3, 'FontSize', textsz)
ylabel('$\gamma$', 'Interpreter','latex', 'FontSize', textsz+5);


subplot(422)
rr = r2;
title(sprintf('Hover %i, MOP %i', Video_bathy(rr).flight, Video(rr).mop))
xmax = Video_bathy(rr).x10';xmin = Video_bathy(rr).x10'; zmax = (Video_bathy(rr).cbathy.z+Video_bathy(rr).cbathy.zerr)';zmin = (Video_bathy(rr).cbathy.z-Video_bathy(rr).cbathy.zerr)';
xmax(isnan(zmax)==1)=[]; zmax(isnan(zmax)==1)=[]; xmin(isnan(zmin)==1)=[]; zmin(isnan(zmin)==1)=[];
patch([xmax fliplr(xmin)], [zmax fliplr(zmin)],[0.6 0.6 0.6], 'EdgeColor', 'none');alpha 0.4
hold on
plot(Video_bathy(rr).x10, Video_bathy(rr).survey.z, 'k', 'LineWidth', 4);
p(1)=plot(Video_bathy(rr).x10, Video_bathy(rr).cbathy.z, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
p(2)=plot(Video_bathy(rr).x10, Video_bathy(rr).cbathy.cbathy_hErr, 'r--', 'LineWidth', 3);
p(3)=plot(Video_bathy(rr).x10, Video_bathy(rr).cbathy.cbathy_gamma, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
p(4)=plot(Video_bathy(rr).x10, Video_bathy(rr).composite.cbathy_nlin,'--','Color', [0 0.4470 0.7410], 'LineWidth', 2);
xlabel('cBathy')
ylabel('Elevation (NAVD88 m)')
set(gca, 'FontSize', textsz)
grid on
xlim([0 400])
ylim([-8 2])
set(gca, 'box', 'on')

subplot(424)
p(1)=plot(Video_bathy(rr).x10, Video_bathy(rr).survey.z, 'k', 'LineWidth', 4);
hold on
p(2)=plot(Video_bathy(rr).x5, Video_bathy(rr).h_avg.lin,  'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
p(3)=plot(Video_bathy(rr).x5, Video_bathy(rr).h_avg.nlin, '--','Color', [0 0.4470 0.7410], 'LineWidth', 2);
p(4)=plot(Video_bathy(rr).x5, Video_bathy(rr).h_avg.bp, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 3);
xlabel('Crest-tracking')
set(gca, 'FontSize', textsz)
ylabel('Elevation (NAVD88 m)')
grid on
xlim([0 400])
ylim([-8 2])

subplot(426)
p(1)=plot(Video_bathy(rr).x10, Video_bathy(rr).survey.z, 'k', 'LineWidth', 4);
hold on
p(2)=plot(Video_bathy(rr).x10, Video_bathy(rr).composite.cbathy_hErr, 'r--', 'LineWidth', 2);
p(2)=plot(Video_bathy(rr).x10, Video_bathy(rr).composite.cbathy_gamma,  'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
p(3)=plot(Video_bathy(rr).x10, Video_bathy(rr).composite.cbathy_nlin, '--','Color', [0 0.4470 0.7410], 'LineWidth', 2);
p(4)=plot(Video_bathy(rr).x10, Video_bathy(rr).composite.cbathyCT, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 3);
xlabel('Composite')
set(gca, 'FontSize', textsz)
ylabel('Elevation (NAVD88 m)')
grid on
xlim([0 400])
ylim([-8 2])

ax = subplot(428);
xticks([]); yticks([])
if Video_bathy(rr).flight < 10
    [fi,fcmap] = imread(fullfile(data_dir, 'timestacks','data', [char(string(Video_bathy(rr).date)) '_' char(Video_bathy(rr).location) '_0' char(string(Video_bathy(rr).flight)) '_' char(string(Video_bathy(rr).mop)) '.png']));
else
    [fi,fcmap] = imread(fullfile(data_dir, 'timestacks','data', [char(string(Video_bathy(rr).date)) '_' char(Video_bathy(rr).location) '_' char(string(Video_bathy(rr).flight)) '_' char(string(Video_bathy(rr).mop)) '.png']));
end 

ax2 = axes('Position',[ax.Position(1) ax.Position(2) ax.Position(3) ax.Position(4)/2]);
image(ax2,Video_bathy(rr).x10, 1,imrotate(cast(fi, 'uint8'),-90))
set(ax2,'ytick', [])
xlabel('Cross-shore Distance (m)', 'FontSize', textsz-4);
xlim([0 400])
   
ax3 = axes('Position',[ax.Position(1) ax.Position(2)+ax.Position(4)/2 ax.Position(3) ax.Position(4)/2]);
plot(ax3, Video_bathy(rr).x10, Video_bathy(rr).gamma_mean,'Color', [0.4660 0.6740 0.1880], 'LineWidth', 3)
hold on
plot(ax3, [100 100], [0 0.5], 'Color', [0.8 0.8 0.8])
plot(ax3, [200 200], [0 0.5], 'Color', [0.8 0.8 0.8])
plot(ax3, [300 300], [0 0.5], 'Color', [0.8 0.8 0.8])
plot(ax3, [400 400], [0 0.5], 'Color', [0.8 0.8 0.8])
ylim(ax3, [0 0.42])
set(ax3,'xtick', [])
xlim(ax3, [0 400])
set(ax2, 'FontSize', textsz)
set(ax3, 'FontSize', textsz)
ylabel('$\gamma$', 'Interpreter','latex', 'FontSize', textsz+5);

