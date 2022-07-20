%% Generate an example bathymetry plot to visualize all estimated bathymetries
% Pick two example hovers
r1 = 4
r2 = 20

figure(1);clf
subplot(421)
rr = r1;
title(sprintf('Hover %i, MOP %i', Video(rr).flight, Video(rr).mop))
xmax = Video(rr).x10';xmin = Video(rr).x10'; zmax = (Video(rr).cbathy.z+Video(rr).cbathy.zerr)';zmin = (Video(rr).cbathy.z-Video(rr).cbathy.zerr)';
xmax(isnan(zmax)==1)=[]; zmax(isnan(zmax)==1)=[]; xmin(isnan(zmin)==1)=[]; zmin(isnan(zmin)==1)=[];
patch([xmax fliplr(xmin)], [zmax fliplr(zmin)],[0.6 0.6 0.6], 'EdgeColor', 'none');alpha 0.4
hold on
p(1)=plot(Video(rr).x10, Video(rr).survey.z, 'k', 'LineWidth', 4);
p(2)=plot(Video(rr).x10, Video(rr).cbathy.z, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
p(3)=plot(Video(rr).x10, Video(rr).cbathy.cbathy_hErr, 'r--', 'LineWidth', 3);
p(4)=plot(Video(rr).x10, Video(rr).cbathy.cbathy_gamma, 'Color', colors(2,:), 'LineWidth', 2);
p(5)=plot(Video(rr).x10, Video(rr).composite.cbathy_nlin,'--','Color', colors(1,:), 'LineWidth', 2);
lg = legend(p, 'Survey','Default', 'hErr < 0.5', 'Breaking', 'Nonlinear', 'Location', 'southwest', 'Box', 'off', 'FontSize', textsz-2);
xlabel('cBathy')
ylabel('Elevation (NAVD88 m)')
set(gca, 'FontSize', textsz)
grid on
xlim([0 400])
ylim([-8 2])
set(gca, 'box', 'on')

subplot(423)
p(1)=plot(Video(rr).x10, Video(rr).survey.z, 'k', 'LineWidth', 4);
hold on
p(2)=plot(Video(rr).x10, Video(rr).h_avg.lin,  'Color', colors(2,:), 'LineWidth', 2);
p(3)=plot(Video(rr).x10, Video(rr).h_avg.nlin, '--','Color', colors(1,:), 'LineWidth', 2);
p(4)=plot(Video(rr).x10, Video(rr).h_avg.bp, 'Color', colors(5,:), 'LineWidth', 3);
lg = legend(p, 'Survey', 'Linear', 'Nonlinear', 'BP', 'Location', 'northeast', 'Box', 'off', 'FontSize', textsz-2, 'Location', 'southwest')
xlabel('Crest-tracking')
set(gca, 'FontSize', textsz)
ylabel('Elevation (NAVD88 m)')
grid on
xlim([0 400])
ylim([-8 2])

subplot(425)
p(1)=plot(Video(rr).x10, Video(rr).survey.z, 'k', 'LineWidth', 4);
hold on
p(2)=plot(Video(rr).x10, Video(rr).composite.cbathy_hErr, 'r--', 'LineWidth', 2);
p(3)=plot(Video(rr).x10, Video(rr).composite.cbathy_gamma,  'Color', colors(2,:), 'LineWidth', 2);
p(4)=plot(Video(rr).x10, Video(rr).composite.cbathy_nlin, '--','Color', colors(1,:), 'LineWidth', 2);
p(5)=plot(Video(rr).x10, Video(rr).composite.cbathyCT, 'Color', colors(5,:), 'LineWidth', 3);
lg = legend(p, 'Survey', 'cB hErr < 0.5', 'cB Breaking', 'cB Nonlinear', 'cBathyCT', 'Location', 'southwest', 'Box', 'off', 'FontSize', textsz-2);
xlabel('Composite')
set(gca, 'FontSize', textsz)
ylabel('Elevation (NAVD88 m)')
grid on
xlim([0 400])
ylim([-8 2])

ax = subplot(427);
xticks([]); yticks([])
if Video(rr).flight < 10
    [fi,fcmap] = imread(fullfile(data_dir, 'timestacks','data', [char(string(Video(rr).date)) '_' char(Video(rr).location) '_0' char(string(Video(rr).flight)) '_' char(string(Video(rr).mop)) '.png']));
else
    [fi,fcmap] = imread(fullfile(data_dir, 'timestacks','data', [char(string(Video(rr).date)) '_' char(Video(rr).location) '_' char(string(Video(rr).flight)) '_' char(string(Video(rr).mop)) '.png']));
end 

ax2 = axes('Position',[ax.Position(1) ax.Position(2) ax.Position(3) ax.Position(4)/2]);
image(ax2,Video(rr).x10, 1,imrotate(cast(fi, 'uint8'),-90))
set(ax2,'ytick', [])
xlabel('Cross-shore Distance (m)', 'FontSize', textsz-4);
xlim([0 400])
   
ax3 = axes('Position',[ax.Position(1) ax.Position(2)+ax.Position(4)/2 ax.Position(3) ax.Position(4)/2]);
aa=find(~isnan(Video(rr).h_avg.bp));
plot(ax3, Video(rr).x10(aa(1):aa(end)), Video(rr).gamma_mean(aa(1):aa(end)),'Color', colors(5,:), 'LineWidth', 3)
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
title(sprintf('Hover %i, MOP %i', Video(rr).flight, Video(rr).mop))
xmax = Video(rr).x10';xmin = Video(rr).x10'; zmax = (Video(rr).cbathy.z+Video(rr).cbathy.zerr)';zmin = (Video(rr).cbathy.z-Video(rr).cbathy.zerr)';
xmax(isnan(zmax)==1)=[]; zmax(isnan(zmax)==1)=[]; xmin(isnan(zmin)==1)=[]; zmin(isnan(zmin)==1)=[];
patch([xmax fliplr(xmin)], [zmax fliplr(zmin)],[0.6 0.6 0.6], 'EdgeColor', 'none');alpha 0.4
hold on
plot(Video(rr).x10, Video(rr).survey.z, 'k', 'LineWidth', 4);
p(1)=plot(Video(rr).x10, Video(rr).cbathy.z, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
p(2)=plot(Video(rr).x10, Video(rr).cbathy.cbathy_hErr, 'r--', 'LineWidth', 3);
p(3)=plot(Video(rr).x10, Video(rr).cbathy.cbathy_gamma, 'Color', colors(2,:), 'LineWidth', 2);
p(4)=plot(Video(rr).x10, Video(rr).composite.cbathy_nlin,'--','Color', colors(1,:), 'LineWidth', 2);
xlabel('cBathy')
ylabel('Elevation (NAVD88 m)')
set(gca, 'FontSize', textsz)
grid on
xlim([0 400])
ylim([-8 2])
set(gca, 'box', 'on')

subplot(424)
p(1)=plot(Video(rr).x10, Video(rr).survey.z, 'k', 'LineWidth', 4);
hold on
p(2)=plot(Video(rr).x10, Video(rr).h_avg.lin,  'Color', colors(2,:), 'LineWidth', 2);
p(3)=plot(Video(rr).x10, Video(rr).h_avg.nlin, '--','Color', colors(1,:), 'LineWidth', 2);
p(4)=plot(Video(rr).x10, Video(rr).h_avg.bp, 'Color', colors(5,:), 'LineWidth', 3);
xlabel('Crest-tracking')
set(gca, 'FontSize', textsz)
ylabel('Elevation (NAVD88 m)')
grid on
xlim([0 400])
ylim([-8 2])

subplot(426)
p(1)=plot(Video(rr).x10, Video(rr).survey.z, 'k', 'LineWidth', 4);
hold on
p(2)=plot(Video(rr).x10, Video(rr).composite.cbathy_hErr, 'r--', 'LineWidth', 2);
p(2)=plot(Video(rr).x10, Video(rr).composite.cbathy_gamma,  'Color', colors(2,:), 'LineWidth', 2);
p(3)=plot(Video(rr).x10, Video(rr).composite.cbathy_nlin, '--','Color', colors(1,:), 'LineWidth', 2);
p(4)=plot(Video(rr).x10, Video(rr).composite.cbathyCT, 'Color', colors(5,:), 'LineWidth', 3);
xlabel('Composite')
set(gca, 'FontSize', textsz)
ylabel('Elevation (NAVD88 m)')
grid on
xlim([0 400])
ylim([-8 2])

ax = subplot(428);
xticks([]); yticks([])
if Video(rr).flight < 10
    [fi,fcmap] = imread(fullfile(data_dir, 'timestacks','data', [char(string(Video(rr).date)) '_' char(Video(rr).location) '_0' char(string(Video(rr).flight)) '_' char(string(Video(rr).mop)) '.png']));
else
    [fi,fcmap] = imread(fullfile(data_dir, 'timestacks','data', [char(string(Video(rr).date)) '_' char(Video(rr).location) '_' char(string(Video(rr).flight)) '_' char(string(Video(rr).mop)) '.png']));
end 

ax2 = axes('Position',[ax.Position(1) ax.Position(2) ax.Position(3) ax.Position(4)/2]);
image(ax2,Video(rr).x10, 1,imrotate(cast(fi, 'uint8'),-90))
set(ax2,'ytick', [])
xlabel('Cross-shore Distance (m)', 'FontSize', textsz-4);
xlim([0 400])
   
ax3 = axes('Position',[ax.Position(1) ax.Position(2)+ax.Position(4)/2 ax.Position(3) ax.Position(4)/2]);
aa=find(~isnan(Video(rr).h_avg.bp));
plot(ax3, Video(rr).x10(aa(1):aa(end)), Video(rr).gamma_mean(aa(1):aa(end)),'Color', colors(5,:), 'LineWidth', 3)
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

