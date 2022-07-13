%% Lange et al. Remote Sensing Paper Plots
%% Load variables

clear all
%load('/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/Video.mat')
load('/Users/athinalange/Desktop/Video.mat')

load '/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/colors.mat'
load '/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/survey.mat'
clear bp tide

for ii = length(Video):-1:1
    if isempty(Video(ii).xshift)
        Video(ii)=[];
    end
end


for rr = 1:length(Video)
    Video(rr).lims=round(Video(rr).lims);
end
date = [Video.date];
location_num = [Video.location_num];
mop = [Video.mop];
flight = [Video.flight];
tide = [Video.tide];
id_Torrey = find(location_num == 1);
id_Cardiff = find(location_num == 2);
id_SIO = find(location_num == 3);
id_20200707 = find(date == 20200707);
id_20210709 = find(date == 20210709);
id_20210712 = find(date == 20210712);
id_20211026 = find(date == 20211026);
id_20211102 = find(date == 20211102);
id_20211202 = find(date == 20211202);
id_20211215 = find(date == 20211215);
id_wholemop = find(rem(mop,1)==0);
id_TP_mop = id_wholemop(find([Video(id_wholemop).location_num] == 1));

date_opt = ['id_20200707';'id_20210709'; 'id_20210712'; 'id_20211026'; 'id_20211102'; 'id_20211202'; 'id_20211215'];

for rr = 1:length(Video)
    Video(rr).Error.std_err = nanstd(Video(rr).h_preavg.bp,[],2)./sqrt(floor(size(Video(rr).h_preavg.bp,2)/10));

    % Extract RMSE and Bias
    RMSE.insz_lin(rr) = Video(rr).Error.RMSE_insz.lin;
    RMSE.insz_nlin(rr) = Video(rr).Error.RMSE_insz.nlin;
    RMSE.insz_bp(rr) = Video(rr).Error.RMSE_insz.bp;
    RMSE.break_lin(rr) = Video(rr).Error.RMSE_break.lin;
    RMSE.break_nlin(rr) = Video(rr).Error.RMSE_break.nlin;
    RMSE.break_bp(rr) = Video(rr).Error.RMSE_break.bp;
    RMSE.sz_lin(rr) = Video(rr).Error.RMSE_sz.lin;
    RMSE.sz_nlin(rr) = Video(rr).Error.RMSE_sz.nlin;
    RMSE.sz_bp(rr) = Video(rr).Error.RMSE_sz.bp;

    Bias.insz_lin(rr) = Video(rr).Error.Bias_insz.lin;
    Bias.insz_nlin(rr) = Video(rr).Error.Bias_insz.nlin;
    Bias.insz_bp(rr) = Video(rr).Error.Bias_insz.bp;
    Bias.break_lin(rr) = Video(rr).Error.Bias_break.lin;
    Bias.break_nlin(rr) = Video(rr).Error.Bias_break.nlin;
    Bias.break_bp(rr) = Video(rr).Error.Bias_break.bp;
    Bias.sz_lin(rr) = Video(rr).Error.Bias_sz.lin;
    Bias.sz_nlin(rr) = Video(rr).Error.Bias_sz.nlin;
    Bias.sz_bp(rr) = Video(rr).Error.Bias_sz.bp;
    
    if isfield(Video(rr).Error.RMSE_insz,'cb')
    
        RMSE.insz_cb(rr) = Video(rr).Error.RMSE_insz.cb;
        RMSE.insz_cb_hErr(rr) = Video(rr).Error.RMSE_insz.cb_hErr;
        RMSE.insz_cb_gamma(rr) = Video(rr).Error.RMSE_insz.cb_gamma;
        RMSE.insz_comp_hErr(rr) = Video(rr).Error.RMSE_insz.comp_hErr;
        RMSE.insz_comp_gamma(rr) = Video(rr).Error.RMSE_insz.comp_gamma;
        RMSE.insz_comp_nlin(rr) = Video(rr).Error.RMSE_insz.comp_nlin; 
        RMSE.insz_comp_bp(rr) = Video(rr).Error.RMSE_insz.comp_bp; 
    
        RMSE.break_cb(rr) = Video(rr).Error.RMSE_break.cb;
        RMSE.break_cb_hErr(rr) = Video(rr).Error.RMSE_break.cb_hErr;
        RMSE.break_cb_gamma(rr) = Video(rr).Error.RMSE_break.cb_gamma;
        RMSE.break_comp_hErr(rr) = Video(rr).Error.RMSE_break.comp_hErr;
        RMSE.break_comp_gamma(rr) = Video(rr).Error.RMSE_break.comp_gamma;
        RMSE.break_comp_nlin(rr) = Video(rr).Error.RMSE_break.comp_nlin; 
        RMSE.break_comp_bp(rr) = Video(rr).Error.RMSE_break.comp_bp; 
        
        RMSE.sz_cb(rr) = Video(rr).Error.RMSE_sz.cb;
        RMSE.sz_cb_hErr(rr) = Video(rr).Error.RMSE_sz.cb_hErr;
        RMSE.sz_cb_gamma(rr) = Video(rr).Error.RMSE_sz.cb_gamma;
        RMSE.sz_comp_hErr(rr) = Video(rr).Error.RMSE_sz.comp_hErr;
        RMSE.sz_comp_gamma(rr) = Video(rr).Error.RMSE_sz.comp_gamma;
        RMSE.sz_comp_nlin(rr) = Video(rr).Error.RMSE_sz.comp_nlin; 
        RMSE.sz_comp_bp(rr) = Video(rr).Error.RMSE_sz.comp_bp;
    
        RMSE.full_comp_hErr(rr) = Video(rr).Error.RMSE_full.comp_hErr;
        RMSE.full_comp_gamma(rr) = Video(rr).Error.RMSE_full.comp_gamma;
        RMSE.full_comp_nlin(rr) = Video(rr).Error.RMSE_full.comp_nlin; 
        RMSE.full_comp_bp(rr) = Video(rr).Error.RMSE_full.comp_bp;
        
        RMSE.offshore_comp_hErr(rr) = Video(rr).Error.RMSE_offshore.comp_hErr;
        RMSE.offshore_comp_gamma(rr) = Video(rr).Error.RMSE_offshore.comp_gamma;
        RMSE.offshore_comp_nlin(rr) = Video(rr).Error.RMSE_offshore.comp_nlin; 
        RMSE.offshore_comp_bp(rr) = Video(rr).Error.RMSE_offshore.comp_bp;


        Bias.insz_cb(rr) = Video(rr).Error.Bias_insz.cb;
        Bias.insz_cb_hErr(rr) = Video(rr).Error.Bias_insz.cb_hErr;
        Bias.insz_cb_gamma(rr) = Video(rr).Error.Bias_insz.cb_gamma;
        Bias.insz_comp_hErr(rr) = Video(rr).Error.Bias_insz.comp_hErr;
        Bias.insz_comp_gamma(rr) = Video(rr).Error.Bias_insz.comp_gamma;
        Bias.insz_comp_nlin(rr) = Video(rr).Error.Bias_insz.comp_nlin; 
        Bias.insz_comp_bp(rr) = Video(rr).Error.Bias_insz.comp_bp; 
    
        Bias.break_cb(rr) = Video(rr).Error.Bias_break.cb;
        Bias.break_cb_hErr(rr) = Video(rr).Error.Bias_break.cb_hErr;
        Bias.break_cb_gamma(rr) = Video(rr).Error.Bias_break.cb_gamma;
        Bias.break_comp_hErr(rr) = Video(rr).Error.Bias_break.comp_hErr;
        Bias.break_comp_gamma(rr) = Video(rr).Error.Bias_break.comp_gamma;
        Bias.break_comp_nlin(rr) = Video(rr).Error.Bias_break.comp_nlin; 
        Bias.break_comp_bp(rr) = Video(rr).Error.Bias_break.comp_bp; 
        
        Bias.sz_cb(rr) = Video(rr).Error.Bias_sz.cb;
        Bias.sz_cb_hErr(rr) = Video(rr).Error.Bias_sz.cb_hErr;
        Bias.sz_cb_gamma(rr) = Video(rr).Error.Bias_sz.cb_gamma;
        Bias.sz_comp_hErr(rr) = Video(rr).Error.Bias_sz.comp_hErr;
        Bias.sz_comp_gamma(rr) = Video(rr).Error.Bias_sz.comp_gamma;
        Bias.sz_comp_nlin(rr) = Video(rr).Error.Bias_sz.comp_nlin; 
        Bias.sz_comp_bp(rr) = Video(rr).Error.Bias_sz.comp_bp;
    
        Bias.full_comp_hErr(rr) = Video(rr).Error.Bias_full.comp_hErr;
        Bias.full_comp_gamma(rr) = Video(rr).Error.Bias_full.comp_gamma;
        Bias.full_comp_nlin(rr) = Video(rr).Error.Bias_full.comp_nlin; 
        Bias.full_comp_bp(rr) = Video(rr).Error.Bias_full.comp_bp;
        
        Bias.offshore_comp_hErr(rr) = Video(rr).Error.Bias_offshore.comp_hErr;
        Bias.offshore_comp_gamma(rr) = Video(rr).Error.Bias_offshore.comp_gamma;
        Bias.offshore_comp_nlin(rr) = Video(rr).Error.Bias_offshore.comp_nlin; 
        Bias.offshore_comp_bp(rr) = Video(rr).Error.Bias_offshore.comp_bp;
    else
        RMSE.insz_cb(rr) = NaN;
        RMSE.insz_cb_hErr(rr) = NaN;
        RMSE.insz_cb_gamma(rr) = NaN;
        RMSE.insz_comp_hErr(rr) = NaN;
        RMSE.insz_comp_gamma(rr) = NaN;
        RMSE.insz_comp_nlin(rr) = NaN; 
        RMSE.insz_comp_bp(rr) = NaN; 
    
        RMSE.break_cb(rr) = NaN;
        RMSE.break_cb_hErr(rr) = NaN;
        RMSE.break_cb_gamma(rr) = NaN;
        RMSE.break_comp_hErr(rr) = NaN;
        RMSE.break_comp_gamma(rr) = NaN;
        RMSE.break_comp_nlin(rr) = NaN; 
        RMSE.break_comp_bp(rr) = NaN; 
        
        RMSE.sz_cb(rr) = NaN;
        RMSE.sz_cb_hErr(rr) = NaN;
        RMSE.sz_cb_gamma(rr) = NaN;
        RMSE.sz_comp_hErr(rr) = NaN;
        RMSE.sz_comp_gamma(rr) = NaN;
        RMSE.sz_comp_nlin(rr) = NaN; 
        RMSE.sz_comp_bp(rr) = NaN;
    
        RMSE.full_comp_hErr(rr) = NaN;
        RMSE.full_comp_gamma(rr) = NaN;
        RMSE.full_comp_nlin(rr) = NaN; 
        RMSE.full_comp_bp(rr) = NaN;
        
        RMSE.offshore_comp_hErr(rr) = NaN;
        RMSE.offshore_comp_gamma(rr) = NaN;
        RMSE.offshore_comp_nlin(rr) = NaN; 
        RMSE.offshore_comp_bp(rr) = NaN;


        Bias.insz_cb(rr) = NaN;
        Bias.insz_cb_hErr(rr) = NaN;
        Bias.insz_cb_gamma(rr) = NaN;
        Bias.insz_comp_hErr(rr) = NaN;
        Bias.insz_comp_gamma(rr) = NaN;
        Bias.insz_comp_nlin(rr) = NaN; 
        Bias.insz_comp_bp(rr) = NaN; 
    
        Bias.break_cb(rr) = NaN;
        Bias.break_cb_hErr(rr) = NaN;
        Bias.break_cb_gamma(rr) = NaN;
        Bias.break_comp_hErr(rr) = NaN;
        Bias.break_comp_gamma(rr) = NaN;
        Bias.break_comp_nlin(rr) = NaN; 
        Bias.break_comp_bp(rr) = NaN; 
        
        Bias.sz_cb(rr) = NaN;
        Bias.sz_cb_hErr(rr) = NaN;
        Bias.sz_cb_gamma(rr) = NaN;
        Bias.sz_comp_hErr(rr) = NaN;
        Bias.sz_comp_gamma(rr) = NaN;
        Bias.sz_comp_nlin(rr) = NaN; 
        Bias.sz_comp_bp(rr) = NaN;
    
        Bias.full_comp_hErr(rr) = NaN;
        Bias.full_comp_gamma(rr) = NaN;
        Bias.full_comp_nlin(rr) = NaN; 
        Bias.full_comp_bp(rr) = NaN;
        
        Bias.offshore_comp_hErr(rr) = NaN;
        Bias.offshore_comp_gamma(rr) = NaN;
        Bias.offshore_comp_nlin(rr) = NaN; 
        Bias.offshore_comp_bp(rr) = NaN;

    end
    
    Swidth(rr) = (Video(rr).lims(3)-Video(rr).lims(2))/10;

end

textsz = 14
r1=find([Video.date]==20211215 & [Video.flight]==1 & [Video.mop]== 583)
r2=find([Video.date]==20211202 & [Video.flight]==1 & [Video.mop]== 582)


%clear id aa rr jj mm ff date_opt x10 x y z zerr ids flight_opt
clearvars -except Video date flight location_num mop tide id_* RMSE* colors Swidth textsz r1 r2 survey Bias*


%% Fig 1 - cBathy vs Timestack Domain
[hFig, ax] = makeFig(7,7,2,2)
rr=r1
% Load timestack
if Video(rr).flight < 10
    gt = imread(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/timestacks/processed/' char(string(Video(rr).date)) '_' char(Video(rr).location) '_0' char(string(Video(rr).flight)) '_' char(string(Video(rr).mop)) '_prediction.jpg']);
    fi = imread(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/timestacks/data/' char(string(Video(rr).date)) '_' char(Video(rr).location) '_0' char(string(Video(rr).flight)) '_' char(string(Video(rr).mop)) '.png']);
else
    gt = imread(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/timestacks/processed/' char(string(Video(rr).date)) '_' char(Video(rr).location) '_' char(string(Video(rr).flight)) '_' char(string(Video(rr).mop)) '_prediction.jpg']);
    fi = imread(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/timestacks/data/' char(string(Video(rr).date)) '_' char(Video(rr).location) '_' char(string(Video(rr).flight)) '_' char(string(Video(rr).mop)) '.png']);
end                          
if size(gt,3) ~= 1
    gt = rgb2gray(gt);
end

% Average line 
gt_clean = bwskel(imbinarize(imcomplement(gt)),'MinBranchLength',20);
gt_gray = ones(size(gt))*255;
gt_gray(gt_clean==1)=0;

aa=gt_gray';aa(aa==255)=NaN;
[X,T]=meshgrid([1:size(aa,2)],[1:size(aa,1)]);

if Video(rr).flight < 10
    cB = imread(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/' char(string(Video(rr).date)) '_' char(Video(rr).location) '/0' char(string(Video(rr).flight)) '/' char(string(Video(rr).date)) '_' char(Video(rr).location) '_0' char(string(Video(rr).flight)) '_Local.png']);
    localX = load(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/' char(string(Video(rr).date)) '_' char(Video(rr).location) '/0' char(string(Video(rr).flight)) '/Processed_data/GRID_' char(string(Video(rr).date)) '_' char(Video(rr).location) '_0' char(string(Video(rr).flight)) '_NAVD88_2Hz.mat'], 'localX');
    localY = load(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/' char(string(Video(rr).date)) '_' char(Video(rr).location) '/0' char(string(Video(rr).flight)) '/Processed_data/GRID_' char(string(Video(rr).date)) '_' char(Video(rr).location) '_0' char(string(Video(rr).flight)) '_NAVD88_2Hz.mat'], 'localY');
else
    cB = imread(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/' char(string(Video(rr).date)) '_' char(Video(rr).location) '/' char(string(Video(rr).flight)) '/' char(string(Video(rr).date)) '_' char(Video(rr).location) '_' char(string(Video(rr).flight)) '_Local.png']);
    localX = load(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/' char(string(Video(rr).date)) '_' char(Video(rr).location) '/' char(string(Video(rr).flight)) '/Processed_data/GRID_' char(string(Video(rr).date)) '_' char(Video(rr).location) '_' char(string(Video(rr).flight)) '_NAVD88_2Hz.mat'], 'localX');
    localY = load(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/' char(string(Video(rr).date)) '_' char(Video(rr).location) '/' char(string(Video(rr).flight)) '/Processed_data/GRID_' char(string(Video(rr).date)) '_' char(Video(rr).location) '_' char(string(Video(rr).flight)) '_NAVD88_2Hz.mat'], 'localY');
end 

%figure(1);clf
ax1 = subplot(221);ax1.Position(1)= ax1.Position(1)+ 0.06;
image(localY.localY(:,1),localX.localX(1,:),imrotate(cB, -90))
hold on
plot([0 0], [-500 0], 'Color', 'w', 'LineWidth', 2)
plot(100*[1 1], [-500 0], 'Color', 'w', 'LineWidth', 2)
plot(-100*[1 1], [-500 0], 'Color', 'w', 'LineWidth', 2)
plot(-200*[1 1], [-500 0], 'Color', 'w', 'LineWidth', 2)
plot(200*[1 1], [-500 0], 'Color', 'w', 'LineWidth', 2)
text(8, -450, '582', 'Color', 'w', 'FontSize', textsz-5)
text(108, -450, '583', 'Color', 'w', 'FontSize', textsz-5)
text(208, -450, '584', 'Color', 'w', 'FontSize', textsz-5)
text(-92, -450, '581', 'Color', 'w', 'FontSize', textsz-5)
text(-192, -450, '580', 'Color', 'w', 'FontSize', textsz-5)
plot(0, -150, 'r.', 'MarkerSize', 20)
plot([100 100], -150+[25 -25], 'Color', colors(2,:), 'LineWidth', 2)
plot([-100 -100], -150+[25 -25], 'Color', colors(2,:), 'LineWidth', 2)
plot([-100 100], -150+[-25 -25], 'Color', colors(2,:), 'LineWidth', 2)
plot([-100 100], -150+[25 25], 'Color', colors(2,:), 'LineWidth', 2)

yticks([-700:100:0])
yticklabels({'700', '600','500', '400', '300', '200', '100', '0'})
ylim([-500 0])
xlim([-400 400])
%set(gca, 'FontSize', textsz)
tl = title('High waves, high tide'); tl.Position(2)=tl.Position(2)+0.1
yl = ylabel('Cross-shore Distance (m)');
p=get(yl, 'Position');p(1) = p(1)-60;
set(yl, 'Position', p)
xl = xlabel('Alongshore Distance (m)');
p=get(xl, 'Position');p(2) = p(2)+25;
set(xl, 'Position', p)
text(-720,-250,'cBathy', 'FontSize', textsz+5, 'Rotation', 90, 'HorizontalAlignment', 'center')
text(-720,400,'Crest-tracking', 'FontSize', textsz+5, 'Rotation', 90, 'HorizontalAlignment', 'center')
text(-380, -470, 'a', 'FontSize', textsz+5, 'Color', 'w')
%
ax1 = subplot(223);ax1.Position(1)= ax1.Position(1)+ 0.06; ax1.Position(2)  = ax1.Position(2) + 0.03;

image(fi)
hold on
scatter3(T(:), X(:),aa(:),8, 'r', 'filled')
scatter(Video(rr).crests.t.*10, Video(rr).crests.x.*10, 8, 'g', 'filled')
xticks([1:601:9000])
xticklabels({'0','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17'})
xlim([600*7 600*9+10])
yticks([0:1000:5000])
yticklabels({'500', '400', '300', '200', '100', '0'})
ylim([0 5001])
%set(gca, 'FontSize', textsz)
%lg = legend('Rejected', 'Accepted');
%lg.Title.String = 'UNet Crests';
xlabel('Time (min)')
yl = ylabel('Cross-shore Distance (m)');
p=get(yl, 'Position');p(1) = p(1)-80;
set(yl, 'Position', p)
text(600*7+20, 300, 'c', 'FontSize', textsz+5, 'Color', 'w')


rr=r2
% Load timestack
if Video(rr).flight < 10
    gt = imread(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/timestacks/processed/' char(string(Video(rr).date)) '_' char(Video(rr).location) '_0' char(string(Video(rr).flight)) '_' char(string(Video(rr).mop)) '_prediction.jpg']);
    fi = imread(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/timestacks/data/' char(string(Video(rr).date)) '_' char(Video(rr).location) '_0' char(string(Video(rr).flight)) '_' char(string(Video(rr).mop)) '.png']);
else
    gt = imread(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/timestacks/processed/' char(string(Video(rr).date)) '_' char(Video(rr).location) '_' char(string(Video(rr).flight)) '_' char(string(Video(rr).mop)) '_prediction.jpg']);
    fi = imread(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/timestacks/data/' char(string(Video(rr).date)) '_' char(Video(rr).location) '_' char(string(Video(rr).flight)) '_' char(string(Video(rr).mop)) '.png']);
end                          
if size(gt,3) ~= 1
    gt = rgb2gray(gt);
end

% Average line 
gt_clean = bwskel(imbinarize(imcomplement(gt)),'MinBranchLength',20);
gt_gray = ones(size(gt))*255;
gt_gray(gt_clean==1)=0;

aa=gt_gray';aa(aa==255)=NaN;
[X,T]=meshgrid([1:size(aa,2)],[1:size(aa,1)]);

if Video(rr).flight < 10
    cB = imread(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/' char(string(Video(rr).date)) '_' char(Video(rr).location) '/0' char(string(Video(rr).flight)) '/' char(string(Video(rr).date)) '_' char(Video(rr).location) '_0' char(string(Video(rr).flight)) '_Local.png']);
    localX = load(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/' char(string(Video(rr).date)) '_' char(Video(rr).location) '/0' char(string(Video(rr).flight)) '/Processed_data/GRID_' char(string(Video(rr).date)) '_' char(Video(rr).location) '_0' char(string(Video(rr).flight)) '_NAVD88_2Hz.mat'], 'localX');
    localY = load(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/' char(string(Video(rr).date)) '_' char(Video(rr).location) '/0' char(string(Video(rr).flight)) '/Processed_data/GRID_' char(string(Video(rr).date)) '_' char(Video(rr).location) '_0' char(string(Video(rr).flight)) '_NAVD88_2Hz.mat'], 'localY');
else
    cB = imread(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/' char(string(Video(rr).date)) '_' char(Video(rr).location) '/' char(string(Video(rr).flight)) '/' char(string(Video(rr).date)) '_' char(Video(rr).location) '_' char(string(Video(rr).flight)) '_Local.png']);
    localX = load(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/' char(string(Video(rr).date)) '_' char(Video(rr).location) '/' char(string(Video(rr).flight)) '/Processed_data/GRID_' char(string(Video(rr).date)) '_' char(Video(rr).location) '_' char(string(Video(rr).flight)) '_NAVD88_2Hz.mat'], 'localX');
    localY = load(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/' char(string(Video(rr).date)) '_' char(Video(rr).location) '/' char(string(Video(rr).flight)) '/Processed_data/GRID_' char(string(Video(rr).date)) '_' char(Video(rr).location) '_' char(string(Video(rr).flight)) '_NAVD88_2Hz.mat'], 'localY');
end 


ax1 = subplot(224);ax1.Position(1)= ax1.Position(1)+ 0.01; ax1.Position(2)  = ax1.Position(2) + 0.03;

image(fi)
hold on
scatter3(T(:), X(:),aa(:),8, 'r', 'filled')
scatter(Video(rr).crests.t.*10, Video(rr).crests.x.*10, 8, 'g', 'filled')
xticks([1:601:9000])
xticklabels({'0','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17'})
xlim([600*7 600*9+10])
xlabel('Time (min)')

yticks([0:1000:5000])
yticklabels({''})
ylim([0 5001])
%set(gca, 'FontSize', textsz)
lg = legend('Rejected', 'Accepted');
lg.Title.String = 'UNet Crests';
text(600*7+20, 300, 'd', 'FontSize', textsz+5, 'Color', 'w')


ax1 = subplot(222);ax1.Position(1)= ax1.Position(1)+ 0.01;

image(localY.localY(:,1),localX.localX(1,:),imrotate(cB, -90))
hold on
plot([0 0], [-500 0], 'Color', 'w', 'LineWidth', 2)
plot(100*[1 1], [-500 0], 'Color', 'w', 'LineWidth', 2)
plot(-100*[1 1], [-500 0], 'Color', 'w', 'LineWidth', 2)
plot(-200*[1 1], [-500 0], 'Color', 'w', 'LineWidth', 2)
plot(200*[1 1], [-500 0], 'Color', 'w', 'LineWidth', 2)
text(8, -450, '582', 'Color', 'w', 'FontSize', textsz-5)
text(108, -450, '583', 'Color', 'w', 'FontSize', textsz-5)
text(208, -450, '584', 'Color', 'w', 'FontSize', textsz-5)
text(-92, -450, '581', 'Color', 'w', 'FontSize', textsz-5)
text(-192, -450, '580', 'Color', 'w', 'FontSize', textsz-5)

plot(0, -150, 'r.', 'MarkerSize', 20)
plot([100 100], -150+[25 -25], 'Color', colors(2,:), 'LineWidth', 2)
plot([-100 -100], -150+[25 -25], 'Color', colors(2,:), 'LineWidth', 2)
plot([-100 100], -150+[-25 -25], 'Color', colors(2,:), 'LineWidth', 2)
plot([-100 100], -150+[25 25], 'Color', colors(2,:), 'LineWidth', 2)
text(-380, -470, 'b', 'FontSize', textsz+5, 'Color', 'w')

yticks([-700:100:0])
yticklabels({''})
ylim([-500 0])
xlim([-400 400])
%set(gca, 'FontSize', textsz)
tl = title('Low waves, low tide'); tl.Position(2)=tl.Position(2)+0.1
xl = xlabel('Alongshore Distance (m)');
p=get(xl, 'Position');p(2) = p(2)+25;
set(xl, 'Position', p)

%set(gcf, 'Position', [1400 500 2025 1750])
%saveas(gcf, ['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/Figures/Crest_v_cBathy_region.pdf'])
set(0,'DefaultFigureColor','remove')
print('/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/Figures/FIG1','-dpdf', '-r300', '-bestfit')
close all
%clearvars -except Video date flight location_num mop tide id_* cbathy* survey* crest_track* RMSE* Bias* colors bp Swidth textsz r1 r2

%% Fig 2 - range of observed bathys
[hFig, ax1] = makeFig(7,7,2,1)

clear p
%figure(1);clf
ax = subplot(211); ax.Position(4) = 0.4; ax.Position(2) = ax.Position(2) - 0.026


ztemp = survey;ztemp(2:3)=[];
for ii = 1:length(ztemp)
     z(1+7*(ii-1):7+7*(ii-1),:)=ztemp(ii).z;
end

p(1) = plot(Video(1).x10, nanmean(z,1), 'k','LineWidth', 4)
hold on
p(5) = plot(Video(1).x10, z(1,:),'Color', [0.7 0.7 0.7],'LineWidth', 4)
plot(Video(1).x10, z,'Color', [0.7 0.7 0.7],'LineWidth', 4)
plot(Video(1).x10, nanmean(z,1), 'k','LineWidth', 4)


c=turbo(7);c(1,:)=[];
p(2) = plot([-100 -10], [-2 2],'LineWidth', 3, 'Color', c(1,:)); % 20200707
hold on
p(3) = plot([-100 -10], [-2 2],'LineWidth', 3, 'Color', c(2,:)); % 20210712
p(4) = plot([-100 -10], [-2 2],'LineWidth', 3, 'Color', c(3,:)); % 20211026
p(6) = plot([-100 -10], [-2 2],'LineWidth', 3, 'Color', c(4,:)); % 20211102
p(7) = plot([-100 -10], [-2 2],'LineWidth', 3, 'Color', c(5,:)); % 20211202
p(8) = plot([-100 -10], [-2 2],'LineWidth', 3, 'Color', c(6,:)); % 20211215
lg = legend(p, 'Mean h(x)', '2020 07/07', '2021 07/12', '2021 10/26', 'All h(x)', '2021 11/02', '2021 12/02', '2021 12/15')
lg.FontSize = textsz-2;
lg.NumColumns = 2;
xlim([0 500])
ylim([-10 5])
grid on
set(gca, 'ytick', [-10:5:5])
title('Torrey Pines', 'FontSize', textsz)
set(gca, 'FontSize', textsz)
grid on
ylabel('Elevation (NAVD88 m)')
xlabel('Cross-shore Distance (m)')
text(30, 4, 'a', 'FontSize', textsz+5, 'Color','k')


% 
ax2 = axes('InnerPosition',[ax.Position(1)+0.048 ax.Position(2)+0.0315 ax.Position(3)*0.4 ax.Position(4)*0.3])
c=turbo(7);c(1,:)=[];
plot(ax2,Video(1).x10, nanmean(z(4:7:end,:),1) - z(4,:),'LineWidth', 3, 'Color', c(1,:)) % 20200707
hold on
plot(ax2,Video(1).x10, nanmean(z(4:7:end,:),1) - z(11,:),'LineWidth', 3, 'Color', c(2,:)) % 20210712
plot(ax2,Video(1).x10, nanmean(z(4:7:end,:),1) - z(18,:),'LineWidth', 3, 'Color', c(3,:)) % 20211026
plot(ax2,Video(1).x10, nanmean(z(4:7:end,:),1) - z(25,:),'LineWidth', 3, 'Color', c(4,:)) % 20211102
plot(ax2,Video(1).x10, nanmean(z(4:7:end,:),1) - z(32,:),'LineWidth', 3, 'Color', c(5,:)) % 20211202
plot(ax2,Video(1).x10, nanmean(z(4:7:end,:),1) - z(39,:),'LineWidth', 3, 'Color', c(6,:)) % 20211215
plot([0 500], [0 0], 'k', 'LineWidth', 3)

title('Mean - survey on 582')
set(ax2, 'FontSize', textsz-4)
grid on
xlim([0 500])
xticks([0:100:500])
xticklabels({'0', '', '', '', '', '500'})


ax = subplot(212); ax.Position(2) = ax.Position(2) - 0.030;

hold on
clear p
for ii = id_Cardiff
   if rem(Video(ii).mop,1)==0
       if Video(ii).mop == 667
            ic=1;
            p(ic)=plot(0, 5, 'Color', colors(ic,:), 'LineWidth', 4);
        elseif Video(ii).mop == 668
            ic=2;
            p(ic)=plot(0, 5, 'Color', colors(ic,:), 'LineWidth', 4);
        elseif Video(ii).mop == 669
            ic=3;
            p(ic)=plot(0, 5, 'Color', colors(ic,:), 'LineWidth', 4);
        elseif Video(ii).mop == 670
            ic=4;
            p(ic)=plot(0, 5, 'Color', colors(ic,:), 'LineWidth', 4);
        end
        plot(Video(ii).x10, Video(ii).survey.z, 'Color', colors(ic,:), 'LineWidth',  4)
   end
end
p(5)=plot(Video(182).x10,Video(182).survey.z, 'Color', colors(5,:), 'LineWidth', 4)
title('Cardiff / SIO', 'FontSize', textsz)
lg = legend(p, '667', '668', '669', '670', 'SIO', 'NumColumns', 2);
lg.Title.String = 'MOP'
set(gca, 'FontSize', textsz)
grid on
ylabel('Elevation (NAVD88 m)')
xlabel('Cross-shore Distance (m)')
text(30, 4, 'b', 'FontSize', textsz+5, 'Color','k')
%set(gcf, 'Position', [0 0 1200 1700])
%saveas(gcf, ['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/Figures/ground_truth_range.png'])

set(0,'DefaultFigureColor','remove')
print('/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/Figures/FIG2.pdf','-dpdf', '-r300', '-bestfit')

%clearvars -except Video date flight location_num mop tide id_* cbathy* survey* crest_track* RMSE* Bias* colors bp Swidth textsz r1 r2
close all

%% Fig 3 - UNET - DONE
%% Fig 4 - Inversion Methods Transects
[hFig, ax1] = makeFig(7,8,3,2)
    
clear p
 
    ax = subplot(321); ax.Position = [ax.Position(1) + 0.04; ax.Position(2) - 0.03; ax.Position(3) + 0.04; ax.Position(4) + 0.06]
    rr=r1
    title('High waves, high tide')
    xmax = Video(rr).x10';xmin = Video(rr).x10'; zmax = (Video(rr).cbathy.z+Video(rr).cbathy.zerr)';zmin = (Video(rr).cbathy.z-Video(rr).cbathy.zerr)';
    xmax(isnan(zmax)==1)=[]; zmax(isnan(zmax)==1)=[]; xmin(isnan(zmin)==1)=[]; zmin(isnan(zmin)==1)=[];
    patch([xmax fliplr(xmin)]+Video(rr).xshift/10, [zmax fliplr(zmin)],[0.6 0.6 0.6], 'EdgeColor', 'none');alpha 0.4
    hold on
    plot(Video(rr).x10+Video(rr).xshift/10, Video(rr).survey.z, 'k', 'LineWidth', 4);
    p(1)=plot(Video(rr).x10+Video(rr).xshift/10, Video(rr).cbathy.z, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    p(2)=plot(Video(rr).x10+Video(rr).xshift/10, Video(rr).cbathy.cbathy_hErr, 'r--', 'LineWidth', 3);
    p(3)=plot(Video(rr).x10+Video(rr).xshift/10, Video(rr).cbathy.cbathy_gamma, 'Color', colors(2,:), 'LineWidth', 2);
    p(4)=plot(Video(rr).x10+Video(rr).xshift/10, Video(rr).composite.cbathy_nlin,'--','Color', colors(1,:), 'LineWidth', 2);
    text(-110, -4, 'cBathy', 'FontSize', textsz + 5, 'HorizontalAlignment', 'center', 'Rotation', 90)
    text(-180, -12.5,'Elevation (NAVD88 m)', 'FontSize', textsz, 'HorizontalAlignment', 'center', 'Rotation', 90)
    set(gca, 'FontSize', textsz)
    grid on
    xlim([0 500])
    ylim([-10 2])
    text(450, 1, 'a', 'Color', 'k', 'FontSize', textsz+5)
    set(ax, 'box', 'on')

    ax = subplot(322); ax.Position = [ax.Position(1) + 0.03; ax.Position(2) - 0.03; ax.Position(3) + 0.04; ax.Position(4) + 0.06]
    rr=r2
    title('Low waves, low tide')
    xmax = Video(rr).x10';xmin = Video(rr).x10'; zmax = (Video(rr).cbathy.z+Video(rr).cbathy.zerr)';zmin = (Video(rr).cbathy.z-Video(rr).cbathy.zerr)';
    xmax(isnan(zmax)==1)=[]; zmax(isnan(zmax)==1)=[]; xmin(isnan(zmin)==1)=[]; zmin(isnan(zmin)==1)=[];
    patch([xmax fliplr(xmin)]+Video(rr).xshift/10, [zmax fliplr(zmin)],[0.6 0.6 0.6], 'EdgeColor', 'none');alpha 0.4
    hold on
    p(1)=plot(Video(rr).x10+Video(rr).xshift/10, Video(rr).survey.z, 'k', 'LineWidth', 4);
    p(2)=plot(Video(rr).x10+Video(rr).xshift/10, Video(rr).cbathy.z, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    p(3)=plot(Video(rr).x10+Video(rr).xshift/10, Video(rr).cbathy.cbathy_hErr, 'r--', 'LineWidth', 3);
    p(4)=plot(Video(rr).x10+Video(rr).xshift/10, Video(rr).cbathy.cbathy_gamma, 'Color', colors(2,:), 'LineWidth', 2);
    p(5)=plot(Video(rr).x10+Video(rr).xshift/10, Video(rr).composite.cbathy_nlin,'--','Color', colors(1,:), 'LineWidth', 2);
    lg = legend(p, 'Survey','Default', 'hErr < 0.5', 'Breaking', 'Nonlinear', 'Location', 'southwest', 'Box', 'off', 'FontSize', textsz-2);
    lg.Position(1) = lg.Position(1) + 0.03; 
    lg.Position(2) = lg.Position(2) + 0.01;
    set(gca, 'FontSize', textsz)
    grid on
    xlim([0 500])
    ylim([-10 2])
    yticks([-10:5:5])
    yticklabels({''})
    text(450, 1, 'b', 'Color', 'k', 'FontSize', textsz+5)
    set(ax, 'box', 'on')

    ax = subplot(323); ax.Position = [ax.Position(1) + 0.04; ax.Position(2) - 0.06; ax.Position(3) + 0.04; ax.Position(4) + 0.06]
    rr=r1
    p(1)=plot(Video(rr).x10+Video(rr).xshift/10, Video(rr).survey.z, 'k', 'LineWidth', 4);
    hold on
    p(2)=plot(Video(rr).x10+Video(rr).xshift/10, Video(rr).h_avg.lin,  'Color', colors(2,:), 'LineWidth', 2);
    p(3)=plot(Video(rr).x10+Video(rr).xshift/10, Video(rr).h_avg.nlin, '--','Color', colors(1,:), 'LineWidth', 2);
    
    p(4)=plot(Video_5(rr).x10+Video_5(rr).xshift/10, Video_5(rr).h_avg.bp, 'Color', [0.8 0.8 0.8], 'LineWidth', 3);
    p(5)=plot(Video(rr).x10+Video(rr).xshift/10, Video(rr).h_avg.bp, 'Color', colors(5,:), 'LineWidth', 3);
    text(-110, -4, 'Crest-tracking', 'FontSize', textsz+5, 'HorizontalAlignment', 'center', 'Rotation', 90)
    set(gca, 'FontSize', textsz)
    grid on
    xlim([0 500])
    ylim([-10 2])
    text(450, 1, 'c', 'Color', 'k', 'FontSize', textsz+5)

    clear p
    ax = subplot(324); ax.Position = [ax.Position(1) + 0.03; ax.Position(2) - 0.06; ax.Position(3) + 0.04; ax.Position(4) + 0.06]
    rr=r2
    plot(Video(rr).x10+Video(rr).xshift/10, Video(rr).survey.z, 'k', 'LineWidth', 4);
    hold on
    p(1)=plot(Video(rr).x10+Video(rr).xshift/10, Video(rr).h_avg.lin, 'Color', colors(2,:), 'LineWidth', 2);
    p(2)=plot(Video(rr).x10+Video(rr).xshift/10, Video(rr).h_avg.nlin, '--','Color', colors(1,:), 'LineWidth', 2);
    p(3)=plot(Video_5(rr).x10+Video_5(rr).xshift/10, Video_5(rr).h_avg.bp, 'Color', [0.8 0.8 0.8], 'LineWidth', 3);
    p(4)=plot(Video(rr).x10+Video(rr).xshift/10, Video(rr).h_avg.bp, 'Color', colors(5,:), 'LineWidth', 3);
    set(gca, 'FontSize', textsz)
    lg = legend(p, 'Linear', 'Nonlinear','BP (5m)', 'BP (25m)', 'Location', 'northeast', 'Box', 'off', 'FontSize', textsz-2, 'Location', 'southwest')
    lg.Position(1) = lg.Position(1) + 0.03; 
    lg.Position(2) = lg.Position(2) + 0.01;
    grid on
    xlim([0 500])
    ylim([-10 2])
    yticks([-10:5:5])
    yticklabels({''})
    text(450, 1, 'd', 'Color', 'k', 'FontSize', textsz+5)

   
    ax = subplot(325); ax.Position = [ax.Position(1) + 0.04; ax.Position(2) - 0.032; ax.Position(3) + 0.04; ax.Position(4)]
    set(ax, 'xtick', [],'ytick', [])
    rr=r1
    if Video(rr).flight < 10
        fi = imread(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/timestacks/data/' char(string(Video(rr).date)) '_' char(Video(rr).location) '_0' char(string(Video(rr).flight)) '_' char(string(Video(rr).mop)) '.png']);
    else
        fi = imread(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/timestacks/data/' char(string(Video(rr).date)) '_' char(Video(rr).location) '_' char(string(Video(rr).flight)) '_' char(string(Video(rr).mop)) '.png']);
    end 
    ax2 = axes('Position',[ax.Position(1) ax.Position(2) ax.Position(3) ax.Position(4)/2]);
    image(ax2,Video(rr).x10, 1,imrotate(cast(fi, 'uint8'),-90))
    set(ax2,'ytick', [])
    xl = xlabel('Cross-shore Distance (m)', 'FontSize', textsz-4);
    p=get(xl, 'Position');p(2) = p(2)+2000;
    set(xl, 'Position', p)
    
   
    ax3 = axes('Position',[ax.Position(1) ax.Position(2)+ax.Position(4)/2 ax.Position(3) ax.Position(4)/2]);
    aa=find(~isnan(Video(rr).h_avg.bp));
    plot(ax3, Video(rr).x10(aa(1):aa(end))+Video(rr).xshift/10, Video(rr).gamma_mean(aa(1):aa(end)),'Color', colors(5,:), 'LineWidth', 3)
    hold on
    plot(ax3, [100 100], [0 0.5], 'Color', [0.8 0.8 0.8])
    plot(ax3, [200 200], [0 0.5], 'Color', [0.8 0.8 0.8])
    plot(ax3, [300 300], [0 0.5], 'Color', [0.8 0.8 0.8])
    plot(ax3, [400 400], [0 0.5], 'Color', [0.8 0.8 0.8])
    ylim(ax3, [0 0.42])
    set(ax3,'xtick', [])
    xlim(ax3, [0 500])
    set(ax2, 'FontSize', textsz)
    set(ax3, 'FontSize', textsz)
    yl = ylabel('$\gamma$', 'Interpreter','latex', 'FontSize', textsz+5);
    p=get(yl, 'Position');p(1) = p(1)-30;
    set(yl, 'Position', p)
    text(450, 0.32, 'e', 'Color', 'k', 'FontSize', textsz+5)

    ax = subplot(326); ax.Position = [ax.Position(1) + 0.03; ax.Position(2) - 0.032; ax.Position(3) + 0.04; ax.Position(4)]
    set(ax, 'xtick', [],'ytick', [])
    rr=r2
    if Video(rr).flight < 10
        fi = imread(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/timestacks/data/' char(string(Video(rr).date)) '_' char(Video(rr).location) '_0' char(string(Video(rr).flight)) '_' char(string(Video(rr).mop)) '.png']);
    else
        fi = imread(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/timestacks/data/' char(string(Video(rr).date)) '_' char(Video(rr).location) '_' char(string(Video(rr).flight)) '_' char(string(Video(rr).mop)) '.png']);
    end 
    ax2 = axes('Position',[ax.Position(1) ax.Position(2) ax.Position(3) ax.Position(4)/2]);
    image(ax2,Video(rr).x10, 1,imrotate(cast(fi, 'uint8'),-90))
    set(ax2,'ytick', [])
    xl = xlabel('Cross-shore Distance (m)', 'FontSize', textsz-4);
    p=get(xl, 'Position');p(2) = p(2)+2000;
    set(xl, 'Position', p)
    
    ax3 = axes('Position',[ax.Position(1) ax.Position(2)+ax.Position(4)/2 ax.Position(3) ax.Position(4)/2]);
    aa=find(~isnan(Video(rr).h_avg.bp));
    plot(ax3, Video(rr).x10(aa(1):aa(end))+Video(rr).xshift/10, Video(rr).gamma_mean(aa(1):aa(end)),'Color', colors(5,:), 'LineWidth', 3)
    hold on
    plot(ax3, [100 100], [0 0.5], 'Color', [0.8 0.8 0.8])
    plot(ax3, [200 200], [0 0.5], 'Color', [0.8 0.8 0.8])
    plot(ax3, [300 300], [0 0.5], 'Color', [0.8 0.8 0.8])
    plot(ax3, [400 400], [0 0.5], 'Color', [0.8 0.8 0.8])
    ylim(ax3, [0 0.42])
    set(ax3,'xtick', [])
    xlim(ax3, [0 500])
    yticklabels(ax3, {''})
    set(ax2, 'FontSize', textsz)
    set(ax3, 'FontSize', textsz)
    text(ax3, 450, 0.32, 'f', 'Color', 'k', 'FontSize', textsz+5)

   
    set(0,'DefaultFigureColor','remove')
    print('/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/Figures/FIG4.pdf','-dpdf', '-r300', '-bestfit')

    close all


%% Fig 5 - Scatter
[hFig, ax1] = makeFig(14,9,2,4)
    

Video_short = Video(id_20200707);
Video_short([Video_short.mop]~=582)=[];
Video_short = Video_short([1,4,6, 9, 11]);
clear p

c = parula(15);c(1,:)=[];c(end,:)=[]; c = c([1,4,7,10,12],:);
    
ax = subplot(241); ax.Position = [0.045; 0.56; 0.22; 0.35]
p(1)=plot(Video_short(1).x10, Video_short(1).survey.z, 'k', 'LineWidth', 4)
hold on
for ii = 1:5
    p(ii+1)=plot(Video_short(ii).x10, Video_short(ii).cbathy.cbathy_hErr, 'Color', c(ii,:), 'LineWidth', 3)
    plot([Video_short(ii).x10(Video_short(ii).lims(3)) Video_short(ii).x10(Video_short(ii).lims(3))], [-4.5 4], '--', 'Color', c(ii,:), 'LineWidth', 2)
end
lgd = legend(p,['Survey', string(round([Video_short(1).tide],2)),...
    string(round([Video_short(2).tide],2)),string(round([Video_short(3).tide],2)),...
    string(round([Video_short(4).tide],2)),string(round([Video_short(5).tide],2))], 'NumColumns', 2, 'Location', 'southwest', 'box', 'off');
lgd.Title.String = 'Tide (m)';

text(-10, -3.5,'Elevation (NAVD88 m)', 'FontSize', textsz, 'HorizontalAlignment', 'center', 'Rotation', 90)  
set(gca, 'FontSize', textsz)
title('cBathy (hErr < 0.5)')
set(gca, 'FontSize', textsz)
ylim([-8 1])
xlim([50 400])
yticks([-8:2:2])
grid on
text(365, 0.5, 'a', 'FontSize', textsz+5)

ax = subplot(242); ax.Position = [0.045+0.24; 0.56; 0.22; 0.35]
plot(Video_short(1).x10, Video_short(1).survey.z, 'k', 'LineWidth', 4);
hold on
for ii = 1:5
    plot(Video_short(ii).x10, Video_short(ii).cbathy.cbathy_gamma, 'Color', c(ii,:), 'LineWidth', 3);
    plot([Video_short(ii).x10(Video_short(ii).lims(3)) Video_short(ii).x10(Video_short(ii).lims(3))], [-10 4], '--', 'Color', c(ii,:), 'LineWidth', 2)
    plot([Video_short(ii).x10(Video_short(ii).lims(2)) Video_short(ii).x10(Video_short(ii).lims(2))], [-10 4], ':', 'Color', c(ii,:), 'LineWidth', 2)
end
title('cBathy (breaking criterion)')
set(gca, 'FontSize', textsz)
ylim([-8 1])
xlim([50 400])
yticks([-8:2:2])
yticklabels({''})
grid on
text(365, 0.5, 'b', 'FontSize', textsz+5)

text(420, -9.4, 'Cross-shore Distance (m)', 'HorizontalAlignment', 'center', 'FontSize', textsz)
text(420, -22, 'Survey (NAVD88 m)', 'HorizontalAlignment', 'center', 'FontSize', textsz)

text(420, 2.5,'Tidal level varies', 'FontSize', textsz+10, 'HorizontalAlignment', 'center')
text(420, -10.7,'Wave height varies', 'FontSize', textsz+10, 'HorizontalAlignment', 'center')


ax = subplot(244); ax.Position = [0.045+0.24+0.24+0.24; 0.56; 0.22; 0.35]
plot(Video_short(1).x10, Video_short(1).survey.z, 'k', 'LineWidth', 4)
hold on
for ii = 1:5
    plot(Video_short(ii).x10, Video_short(ii).composite.bp, 'Color', c(ii,:), 'LineWidth', 3)
    plot([Video_short(ii).x10(Video_short(ii).lims(3)) Video_short(ii).x10(Video_short(ii).lims(3))], [-10 4], '--', 'Color', c(ii,:), 'LineWidth', 2)
end
title('Composite cBathyCT')
set(gca, 'FontSize', textsz)
ylim([-8 1])
xlim([50 400])
yticks([-8:2:2])
yticklabels({''})
grid on
text(365, 0.5, 'd', 'FontSize', textsz+5)

ax = subplot(243); ax.Position = [0.045+0.24+0.24; 0.56; 0.22; 0.35]
plot(Video_short(1).x10, Video_short(1).survey.z, 'k', 'LineWidth', 4)
hold on
for ii = 1:5
    plot(Video_short(ii).x10, Video_short(ii).h_avg.bp, 'Color', c(ii,:), 'LineWidth', 3)
    plot([Video_short(ii).x10(Video_short(ii).lims(3)) Video_short(ii).x10(Video_short(ii).lims(3))], [-10 4], '--', 'Color', c(ii,:), 'LineWidth', 2)
end
title('Crest-tracking BP')
set(gca, 'FontSize', textsz)
ylim([-8 1])
xlim([50 400])
yticks([-8:2:2])
yticklabels({''})
grid on
text(365, 0.5, 'c', 'FontSize', textsz+5)


Video_short = Video(id_TP_mop);
    for rr = length(Video_short):-1:1
        if ~isfield(Video_short(rr).cbathy, 'cbathy_hErr')
            Video_short(rr)=[];
        end
    end
Video_short([Video_short.mop]~=582)=[];
c=turbo(7);c(1,:)=[];c=c(6:-1:1,:)
clear p
ids = ['20211026'; '20211215'; '20211102'; '20200707'; '20211202';'20210712' ]


ax = subplot(245); ax.Position = [0.045; 0.07; 0.22; 0.35]
jj=0;
hold on
for ii = [18, 37, 23, 7, 28, 17]
    jj=jj+1;
%     if jj == 2 | jj == 3
%         plot([Video_short(ii).survey.z(Video_short(ii).lims(3)) Video_short(ii).survey.z(Video_short(ii).lims(3))], [-10 -2.1],'--', 'Color', c(jj,:), 'LineWidth', 2)
%     else
%         plot([Video_short(ii).survey.z(Video_short(ii).lims(3)) Video_short(ii).survey.z(Video_short(ii).lims(3))], [-10 5],'--', 'Color', c(jj,:), 'LineWidth', 2)
%     end
    p(jj)=scatter(Video_short(ii).survey.z, Video_short(ii).cbathy.cbathy_hErr, 20, c(jj,:), 'filled');
end
plot([-8 2], [-8 2], 'k', 'LineWidth', 4)

leg = legend(p, '2.15', '1.92', '1.12','0.86', '0.65', '0.6', 'Location', 'northeast', 'NumColumns', 2, 'box', 'off');
title(leg, 'Hs (m)')
text(2.55, -3,'Estimated (NAVD88 m)', 'FontSize', textsz, 'HorizontalAlignment', 'center', 'Rotation', 90)  
set(gca, 'FontSize', textsz)
set(gca, 'XDir', 'reverse')
ylim([-8 1])
xlim([-8 1])
yticks([-8:2:2])
grid on
text(-7, 0.5, 'e', 'FontSize', textsz+5)
%xlabel('Survey (NAVD88 m)')
set(ax, 'box', 'on')


ax = subplot(246); ax.Position = [0.045+0.24; 0.07; 0.22; 0.35]
hold on
jj=0
for ii = [18, 37, 23, 7, 28, 17]
    jj=jj+1;
%     plot([Video_short(ii).survey.z(Video_short(ii).lims(3)) Video_short(ii).survey.z(Video_short(ii).lims(3))], [-10 5],'--', 'Color', c(jj,:), 'LineWidth', 2)
    scatter(Video_short(ii).survey.z, Video_short(ii).cbathy.cbathy_gamma, 20, c(jj,:), 'filled');
end
plot([-8 2], [-8 2], 'k', 'LineWidth', 4)
set(gca, 'FontSize', textsz)
set(gca, 'XDir', 'reverse')
ylim([-8 1])
xlim([-8 1])
yticks([-8:2:2])
yticklabels({''})
grid on
text(-7, 0.5, 'f', 'FontSize', textsz+5)
%xlabel('Survey (NAVD88 m)')
set(ax, 'box', 'on')

ax = subplot(247); ax.Position = [0.045+0.24+0.24; 0.07; 0.22; 0.35]
hold on
jj=0
for ii = [18, 37, 23, 7, 28, 17]
    jj=jj+1;
%     plot([Video_short(ii).survey.z(Video_short(ii).lims(3)) Video_short(ii).survey.z(Video_short(ii).lims(3))], [-10 5],'--', 'Color', c(jj,:), 'LineWidth', 2)
    scatter(Video_short(ii).survey.z, Video_short(ii).h_avg.bp, 20, c(jj,:), 'filled');
end
plot([-8 2], [-8 2], 'k', 'LineWidth', 4)
set(gca, 'FontSize', textsz)
set(gca, 'XDir', 'reverse')
ylim([-8 1])
xlim([-8 1])
yticks([-8:2:2])
yticklabels({''})
grid on
text(-7, 0.5, 'g', 'FontSize', textsz+5)
%xlabel('Survey (NAVD88 m)')
set(ax, 'box', 'on')

ax = subplot(248); ax.Position = [0.045+0.24+0.24+0.24; 0.07; 0.22; 0.35]
hold on
jj=0
for ii = [18, 37, 23, 7, 28, 17]
    jj=jj+1;
%     plot([Video_short(ii).survey.z(Video_short(ii).lims(3)) Video_short(ii).survey.z(Video_short(ii).lims(3))], [-10 5],'--', 'Color', c(jj,:), 'LineWidth', 2)
    scatter(Video_short(ii).survey.z, Video_short(ii).composite.bp, 20, c(jj,:), 'filled');
end
plot([-8 2], [-8 2], 'k', 'LineWidth', 4)
set(gca, 'FontSize', textsz)
set(gca, 'XDir', 'reverse')
ylim([-8 1])
xlim([-8 1])
yticks([-8:2:2])
yticklabels({''})
grid on
text(-7, 0.5, 'h', 'FontSize', textsz+5)
%xlabel('Survey (NAVD88 m)')
set(ax, 'box', 'on')


set(0,'DefaultFigureColor','remove')
print('/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/Figures/FIG5.pdf','-dpdf', '-r300', '-bestfit')

close all


%% Fig 6 - Error vs Depth
% Video_short = Video(id_Torrey);
% for rr = length(Video_short):-1:1
%     if isempty(Video_short(rr).composite)
%         Video_short(rr)=[];
%     end
% end
% 
% for rr = 1:length(Video_short)
%     obs(rr,:) = Video_short(rr).survey.z;
%     cb_hErr(rr,:) = Video_short(rr).composite.cbathy_hErr;
%     cb_nlin(rr,:) = Video_short(rr).composite.cbathy_nlin;
%     bp(rr,:) = Video_short(rr).composite.bp;
% end
% 
% for rr = 1:length(Video_short)
%     for ii = 1:5001
%         r_hErr(rr,ii) = sqrt(sum((cb_hErr(rr,ii)-obs(rr,ii)).^2)/length(obs(rr,ii)));
%         r_nlin(rr,ii) = sqrt(sum((cb_nlin(rr,ii)-obs(rr,ii)).^2)/length(obs(rr,ii)));
%         r_bp(rr,ii) = sqrt(sum((bp(rr,ii)-obs(rr,ii)).^2)/length(obs(rr,ii)));
%         
%         b_hErr(rr,ii) = sum(cb_hErr(rr,ii) - obs(rr,ii)) / length(obs(rr,ii));
%         b_nlin(rr,ii) = sum(cb_nlin(rr,ii) - obs(rr,ii)) / length(obs(rr,ii));
%         b_bp(rr,ii) = sum(bp(rr,ii) - obs(rr,ii)) / length(obs(rr,ii));
%     end
% end
% save('/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/depth_values.mat', 'obs', 'cb_hErr', 'cb_nlin', 'bp', 'r_hErr', 'r_nlin', 'r_bp', 'b_hErr', 'b_nlin', 'b_bp')
%

load('/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/depth_values.mat')

for ii = 2:-1:-10
    jj=jj+1;
    if ii >= 0
        eval(['id_' char(num2str(ii)) ' = find(obs <= ii & obs > ii-1);']);
    else
        eval(['id_n' char(num2str(abs(ii))) ' = find(obs <= ii & obs > ii-1);']);
    end

end


for ii = 2:-1:-10
    if ii >= 0
        eval(['r_hErr_bin_' char(num2str(ii)) ' = r_hErr(id_' char(num2str(ii)) ');']);
        eval(['r_nlin_bin_' char(num2str(ii)) ' = r_nlin(id_' char(num2str(ii)) ');']);
        eval(['r_bp_bin_' char(num2str(ii)) ' = r_bp(id_' char(num2str(ii)) ');']);
        eval(['b_hErr_bin_' char(num2str(ii)) ' = b_hErr(id_' char(num2str(ii)) ');']);
        eval(['b_nlin_bin_' char(num2str(ii)) ' = b_nlin(id_' char(num2str(ii)) ');']);
        eval(['b_bp_bin_' char(num2str(ii)) ' = b_bp(id_' char(num2str(ii)) ');']);
    else
        eval(['r_hErr_bin_n' char(num2str(abs(ii))) ' = r_hErr(id_n' char(num2str(abs(ii))) ');']);
        eval(['r_nlin_bin_n' char(num2str(abs(ii))) ' = r_nlin(id_n' char(num2str(abs(ii))) ');']);
        eval(['r_bp_bin_n' char(num2str(abs(ii))) ' = r_bp(id_n' char(num2str(abs(ii))) ');']);
        eval(['b_hErr_bin_n' char(num2str(abs(ii))) ' = b_hErr(id_n' char(num2str(abs(ii))) ');']);
        eval(['b_nlin_bin_n' char(num2str(abs(ii))) ' = b_nlin(id_n' char(num2str(abs(ii))) ');']);
        eval(['b_bp_bin_n' char(num2str(abs(ii))) ' = b_bp(id_n' char(num2str(abs(ii))) ');']);
    end
end
jj=0;
for ii = 2:-1:-10
    jj=jj+1;
    if ii >= 0
        eval(['r_hErr_bin(jj) = nanmean(r_hErr_bin_' char(num2str(ii)) ');']);
        eval(['r_nlin_bin(jj) = nanmean(r_nlin_bin_' char(num2str(ii)) ');']);
        eval(['r_bp_bin(jj) = nanmean(r_bp_bin_' char(num2str(ii)) ');']);
        eval(['b_hErr_bin(jj) = nanmean(b_hErr_bin_' char(num2str(ii)) ');']);
        eval(['b_nlin_bin(jj) = nanmean(b_nlin_bin_' char(num2str(ii)) ');']);
        eval(['b_bp_bin(jj) = nanmean(b_bp_bin_' char(num2str(ii)) ');']);

        eval(['r_hErr_bin_std(jj) = nanstd(r_hErr_bin_' char(num2str(ii)) ');']);
        eval(['r_nlin_bin_std(jj) = nanstd(r_nlin_bin_' char(num2str(ii)) ');']);
        eval(['r_bp_bin_std(jj) = nanstd(r_bp_bin_' char(num2str(ii)) ');']);
        eval(['b_hErr_bin_std(jj) = nanstd(b_hErr_bin_' char(num2str(ii)) ');']);
        eval(['b_nlin_bin_std(jj) = nanstd(b_nlin_bin_' char(num2str(ii)) ');']);
        eval(['b_bp_bin_std(jj) = nanstd(b_bp_bin_' char(num2str(ii)) ');']);
        
    else
        eval(['r_hErr_bin(jj) = nanmean(r_hErr_bin_n' char(num2str(abs(ii))) ');']);
        eval(['r_nlin_bin(jj) = nanmean(r_nlin_bin_n' char(num2str(abs(ii))) ');']);
        eval(['r_bp_bin(jj) = nanmean(r_bp_bin_n' char(num2str(abs(ii))) ');']);
        eval(['b_hErr_bin(jj) = nanmean(b_hErr_bin_n' char(num2str(abs(ii))) ');']);
        eval(['b_nlin_bin(jj) = nanmean(b_nlin_bin_n' char(num2str(abs(ii))) ');']);
        eval(['b_bp_bin(jj) = nanmean(b_bp_bin_n' char(num2str(abs(ii))) ');']);

        eval(['r_hErr_bin_std(jj) = nanstd(r_hErr_bin_n' char(num2str(abs(ii))) ');']);
        eval(['r_nlin_bin_std(jj) = nanstd(r_nlin_bin_n' char(num2str(abs(ii))) ');']);
        eval(['r_bp_bin_std(jj) = nanstd(r_bp_bin_n' char(num2str(abs(ii))) ');']);
        eval(['b_hErr_bin_std(jj) = nanstd(b_hErr_bin_n' char(num2str(abs(ii))) ');']);
        eval(['b_nlin_bin_std(jj) = nanstd(b_nlin_bin_n' char(num2str(abs(ii))) ');']);
        eval(['b_bp_bin_std(jj) = nanstd(b_bp_bin_n' char(num2str(abs(ii))) ');']);
    end
end


[hFig, ax1] = makeFig(7,6,1,2)

ax = subplot(211); ax.Position = [ax.Position(1); 0.55; 0.8; 0.4]
errorbar(-[0:-1:-10], b_hErr_bin(3:end), b_hErr_bin_std(3:end), 'r', 'LineWidth', 3)
hold on
errorbar(-[0:-1:-10], b_nlin_bin(3:end), b_nlin_bin_std(3:end), 'Color', colors(1,:), 'LineWidth', 3)
errorbar(-[0:-1:-10], b_bp_bin(3:end), b_bp_bin_std(3:end), 'Color', colors(5,:), 'LineWidth', 3)
plot(-[0:-1:-10], b_hErr_bin(3:end), 'r*-', 'LineWidth', 3, 'MarkerSize', 10)

plot(-[0:-1:-10], b_nlin_bin(3:end), '*-', 'Color', colors(1,:), 'LineWidth', 3, 'MarkerSize', 10)
plot(-[0:-1:-10], b_bp_bin(3:end), '*-', 'Color', colors(5,:), 'LineWidth', 5, 'MarkerSize', 10)
set(gca, 'FontSize', textsz)
plot([-2 10], [0 0], 'k', 'LineWidth', 3)
grid on
lgd = legend('cBathy (hErr < 0.5)', 'cBathy (nonlinear)', 'cBathyCT', 'Location', 'southeast', 'FontSize', textsz);
%lgd.Title.String = 'Composite';
ylabel('$\langle$ Bias $\rangle$ (m)', 'Interpreter', 'latex')
xlim([-0.5 9.5])
xticks([0:1:10])
xticklabels({''})
text(4.5, 2.4, 'Error Dependence on Depth', 'FontSize', textsz+5, 'HorizontalAlignment', 'center')
text(-0.4, 1.7, 'a', 'FontSize', textsz+5)


ax = subplot(212); ax.Position = [ax.Position(1); 0.1; 0.8; 0.4]
errorbar(-[0:-1:-10], r_hErr_bin(3:end), r_hErr_bin_std(3:end), 'r', 'LineWidth', 3)
hold on
errorbar(-[0:-1:-10], r_nlin_bin(3:end), r_nlin_bin_std(3:end), 'Color', colors(1,:), 'LineWidth', 3)
errorbar(-[0:-1:-10], r_bp_bin(3:end), r_bp_bin_std(3:end), 'Color', colors(5,:), 'LineWidth', 3)

plot(-[0:-1:-10], r_hErr_bin(3:end), 'r*-', 'LineWidth', 3, 'MarkerSize', 10)
plot(-[0:-1:-10], r_nlin_bin(3:end), '*-', 'Color', colors(1,:), 'LineWidth', 3, 'MarkerSize', 10)
plot(-[0:-1:-10], r_bp_bin(3:end), '*-', 'Color', colors(5,:), 'LineWidth', 5, 'MarkerSize', 10)
set(gca, 'FontSize', textsz)
grid on
ylabel('$\langle$ RMSE $\rangle$ (m)', 'Interpreter', 'latex')
xlabel('Depth (m)')
xticks([0:1:10])
xticklabels({'0-1', '1-2', '2-3', '3-4', '4-5', '5-6', '6-7', '7-8', '8-9', '9-10'})
xlim([-0.5 9.5])
ylim([0 2])
text(-0.4, 1.8, 'b', 'FontSize', textsz+5)


set(0,'DefaultFigureColor','remove')
print('/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/Figures/FIG6.pdf','-dpdf', '-r300', '-bestfit')

close all


%% Fig 7 - Timestacks
Video_short = Video(id_20200707);
Video_short([Video_short.mop]~=582)=[];
Video_short = Video_short([1,4,6, 9, 11]);
[hFig, ax1] = makeFig(7,5,1,1)
ax1.Position = [ax1.Position(1); ax1.Position(2); 0.85;0.85]
    clear bounds
    c = parula(15);c(1,:)=[];c(end,:)=[]; c = c([1,4,7,10,12],:);
    
    
    for rr = 1:5
        P(rr,:) = InterX([Video_short(rr).x10'; Video_short(rr).survey.z'], [Video_short(rr).x10'; Video_short(rr).tide * ones(1,5001)]);
    end
    for rr = 1
        ['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/timestacks/data/20200707_Torrey_0' char(string(rr)) '_582.png']
        fi = imread(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/timestacks/data/20200707_Torrey_0' char(string(rr)) '_582.png']);
        id=find(Video_short(rr).gamma_mean > 0.02);id = id(end)
        image(fi)
        hold on
        plot([0 10000], 5001-[id(1) id(1)]-250, 'r-', 'LineWidth', 3)
        %plot([0 10000], 5001-[id(1) id(1)],'--', 'Color', colors(2,:), 'LineWidth', 3)
        plot([0 10000], 5001-[round(P(1,1)*10) round(P(1,1)*10)],'Color', c(1,:), 'LineWidth', 3)
        idsize = size(fi,2);
        totsize = idsize;
        plot([totsize totsize], [0 5000], 'g-', 'LineWidth', 2)
        tideid(rr) = string(round(Video_short(rr).tide,2));
        bounds(rr) = totsize;
    end

    ticks(1)=round(totsize/2);

    rr=1;
    for jj = [4,7,10,12]
        rr=rr+1
        if jj < 10
            ['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/timestacks/data/20200707_Torrey_0' char(string(jj)) '_582.png']
            fi = imread(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/timestacks/data/20200707_Torrey_0' char(string(jj)) '_582.png']);
        else
            ['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/timestacks/data/20200707_Torrey_' char(string(jj)) '_582.png']
            fi = imread(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/timestacks/data/20200707_Torrey_' char(string(jj)) '_582.png']) ;
        end 

        id=find(Video_short(rr).gamma_mean > 0.02);id = id(end)
        image([totsize:totsize + size(fi,2)-1], [1:5001],  fi)
        
        hold on
        plot([totsize totsize + size(fi,2)-1], 5001-[id(1) id(1)]-250, 'r-', 'LineWidth', 3)
        %plot([totsize totsize + size(fi,2)-1], 5001-[id(1) id(1)], '--','Color', colors(2,:), 'LineWidth', 3)
        plot([totsize totsize + size(fi,2)-1], 5001-[round(P(rr,1)*10) round(P(rr,1)*10)],'Color', c(rr,:), 'LineWidth', 3)

        idsize = size(fi,2);
        totsize = totsize+idsize;
        plot([totsize totsize], [0 5000], 'g-', 'LineWidth', 2)
        ticks(rr) = totsize - round(idsize/2);
        bounds(rr) = totsize;
        tideid(rr) = string(round(Video_short(rr).tide,2));
    end
    
    xlim([0 totsize])
    yticks([0:1000:5000])
    yticklabels({'500', '400', '300', '200', '100', '0'})
    xticks(ticks)
    xticklabels({'1', '4', '7', '10', '12'});
    ylim([1999 5001])
    set(gca, 'FontSize', textsz)
    xlabel('Hover')
    ylabel('Cross-shore Distance (m)', 'FontSize', textsz)
    title('Torrey Pines - changing tide level')
    
set(0,'DefaultFigureColor','remove')
print('/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/Figures/FIG7.pdf','-dpdf', '-r300', '-bestfit')

close all


%% Fig 9 - Cardiff

[hFig, ax1] = makeFig(6,4,1,2)

Video_Cardiff = Video(id_Cardiff);
Video_Cardiff = Video_Cardiff(find(rem([Video_Cardiff.mop],1)==0));


ax = subplot(121); ax.Position = [0.08; 0.15; 0.39; 0.73]
set(ax, 'box', 'on')

plot(Video_Cardiff(1).x10, Video_Cardiff(1).survey.z, 'k', 'LineWidth', 5)
hold on
plot(Video_Cardiff(1).x10, Video_Cardiff(1).composite.cbathy_hErr, 'Color', 'r', 'LineWidth', 2)
plot(Video_Cardiff(1).x10, Video_Cardiff(1).composite.cbathy_nlin, 'Color', colors(1,:), 'LineWidth', 2)
plot(Video_Cardiff(1).x10, Video_Cardiff(1).composite.bp, 'Color', colors(5,:), 'LineWidth', 2)
text(220, -1.5, '667', 'FontSize', textsz)

offset = 5
plot(Video_Cardiff(1).x10, Video_Cardiff(2).survey.z - offset, 'k', 'LineWidth', 5)
plot(Video_Cardiff(1).x10, Video_Cardiff(2).composite.cbathy_hErr - offset, 'Color', 'r', 'LineWidth', 2)
plot(Video_Cardiff(1).x10, Video_Cardiff(2).composite.cbathy_nlin - offset, 'Color', colors(1,:), 'LineWidth', 2)
plot(Video_Cardiff(1).x10, Video_Cardiff(2).composite.bp - offset, 'Color', colors(5,:), 'LineWidth', 2)
text(220, -7.1, '668', 'FontSize', textsz)

plot(Video_Cardiff(1).x10, Video_Cardiff(3).survey.z - 2*offset, 'k', 'LineWidth', 5)
plot(Video_Cardiff(1).x10, Video_Cardiff(3).composite.cbathy_hErr - 2*offset, 'Color', 'r', 'LineWidth', 2)
plot(Video_Cardiff(1).x10, Video_Cardiff(3).composite.cbathy_nlin - 2*offset, 'Color', colors(1,:), 'LineWidth', 2)
plot(Video_Cardiff(1).x10, Video_Cardiff(3).composite.bp - 2*offset, 'Color', colors(5,:), 'LineWidth', 2)
text(220, -12.2, '669', 'FontSize', textsz)

grid on
xlim([0 350])
ylim([-19 3])
xticks([0:100:400])
yticks([-15:5:5])

text(315, 1.7, 'a', 'FontSize', textsz + 5)

legend('Survey', 'cBathy (hErr < 0.5)', 'cBathy (nonlinear)', 'cBathyCT', 'FontSize', textsz - 6, 'Location', 'southwest', 'box', 'off')
set(gca, 'FontSize', textsz)
title('Cardiff transects', 'FontSize', textsz)
xlabel('Cross-shore Distance (m)', 'FontSize', textsz)
ylabel('Elevation (NAVD88 m)', 'FontSize', textsz)

ax = subplot(122); ax.Position = [0.65; 0.15; 0.34; 0.73]
set(ax, 'box', 'on')
rr=11
fi = imread(['/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/timestacks/data/20210709_Cardiff_04_667.png']);

image(fi)
hold on
scatter(Video_Cardiff(rr).crests.t.*10, Video_Cardiff(rr).crests.x.*10, 2, 'g', 'filled')
xticks([1:601:9000])
xticklabels({'0','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17'})
xlim([600*7 600*8.5+10])
yticks([0:1000:5000])
yticklabels({'500', '400', '300', '200', '100', '0'})
ylim([0 5001])
set(gca, 'FontSize', textsz)
xlabel('Time (min)', 'FontSize',textsz)
yl = ylabel('Cross-shore Distance (m)', 'FontSize', textsz);
p=get(yl, 'Position');p(1) = p(1)-90;
set(yl, 'Position', p)
text(600*8.3+20, 300, 'b', 'FontSize', textsz+5, 'Color', 'w')
title('Cardiff Reef', 'FontSize', textsz)


set(0,'DefaultFigureColor','remove')
print('/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/Figures/FIG9.pdf','-dpdf', '-r300', '-bestfit')

close all
%% Fig 8 - RMSE 1 vs 3 flight
diff_flight_num_avg_v2

%h=histogram(RMSE_single_full);h.BinWidth = 0.25;h7 = h.Values;
%h7=h7/sum(h7); x7 = 0.25*[1:length(h7)];
date=[Video.date];
id_20200707 = find(date == 20200707);
id_20210709 = find(date == 20210709);
id_20210712 = find(date == 20210712);
id_20211026 = find(date == 20211026);
id_20211102 = find(date == 20211102);
id_20211202 = find(date == 20211202);
id_20211215 = find(date == 20211215);
date_opt = ['id_20210712'; 'id_20211202'; 'id_20200707'; 'id_20211102';  'id_20211215';  'id_20211026'];
x=[0:0.25:6.5];
xb = [-5:0.25:5];
    
hcb = zeros(6,length(x)); hnlin = zeros(6,length(x)); hbp = zeros(6,length(x)); bias = zeros(6,length(xb));
hcb3 = zeros(6,length(x)); hnlin3 = zeros(6,length(x)); hbp3 = zeros(6,length(x)); bias3 = zeros(6,length(xb));

    figure(2)
    for rr = 1:size(date_opt,1)
        eval(['id_run = ' date_opt(rr,:) ';'])
     
        h = histogram(RMSE_hErr_full(id_run));h.BinWidth = 0.25;
        hcb(rr,find(h.BinEdges(1)==x):find(h.BinEdges(1)==x)+length(h.BinEdges)-2) = h.Values; clear h
        
        h = histogram(RMSE_nlin_full(id_run));h.BinWidth = 0.25;
        hnlin(rr,find(h.BinEdges(1)==x):find(h.BinEdges(1)==x)+length(h.BinEdges)-2) = h.Values; clear h
        
        h = histogram(RMSE_single_full(id_run));h.BinWidth = 0.25;
        hbp(rr,find(h.BinEdges(1)==x):find(h.BinEdges(1)==x)+length(h.BinEdges)-2) = h.Values; clear h

        h = histogram(Bias_single_full(id_run));h.BinWidth = 0.25;
        bias(rr,find(h.BinEdges(1)==xb):find(h.BinEdges(1)==xb)+length(h.BinEdges)-2) = h.Values; clear h
    end

    
hbp = hbp./nansum(nansum(hbp));
hcb= hcb./nansum(nansum(hcb));
hnlin = hnlin./nansum(nansum(hnlin));
bias = bias./nansum(nansum(bias));

date = [Video_3_bp.date];
id_20200707 = find(date == 20200707);
id_20210709 = find(date == 20210709);
id_20210712 = find(date == 20210712);
id_20211026 = find(date == 20211026);
id_20211102 = find(date == 20211102);
id_20211202 = find(date == 20211202);
id_20211215 = find(date == 20211215);

    figure(2)
    for rr = 1:size(date_opt,1)
        eval(['id_run = ' date_opt(rr,:) ';'])
     
        h = histogram([Video_3_cb_hErr(id_run).RMSE_median_full]);h.BinWidth = 0.25;
        hcb3(rr,find(h.BinEdges(1)==x):find(h.BinEdges(1)==x)+length(h.BinEdges)-2) = h.Values; clear h
        
        h = histogram([Video_3_cb_nlin(id_run).RMSE_median_full]);h.BinWidth = 0.25;
        hnlin3(rr,find(h.BinEdges(1)==x):find(h.BinEdges(1)==x)+length(h.BinEdges)-2) = h.Values; clear h
        
        h = histogram([Video_3_bp(id_run).RMSE_median_full]);h.BinWidth = 0.25;
        hbp3(rr,find(h.BinEdges(1)==x):find(h.BinEdges(1)==x)+length(h.BinEdges)-2) = h.Values; clear h
        
        h = histogram([Video_3_bp(id_run).Bias_median_full]);h.BinWidth = 0.25;
        bias3(rr,find(h.BinEdges(1)==xb):find(h.BinEdges(1)==xb)+length(h.BinEdges)-2) = h.Values; clear h
        
    end

    
hbp3 = hbp3./nansum(nansum(hbp3));
hcb3 = hcb3./nansum(nansum(hcb3));
hnlin3 = hnlin3./nansum(nansum(hnlin3));
bias3 = bias3./nansum(nansum(bias3));



textsz = 14
[hFig, ax1] = makeFig(7,8,4,2)
c=turbo(7);c(1,:)=[];

ax = subplot(421); ax.Position = [0.17; 0.77; 0.4; 0.2]
h=bar(x,hcb,'stacked'); for ii = 1:6;h(ii).FaceColor = c(ii,:);end
tl = title('1 Flight', 'FontSize', textsz);
lgd=legend('0.6m (80m)', '0.65m (112m)', '0.86m (125m)', '1.12m (170m)', '1.92m (225m)', '2.15m (240m)', 'NumColumns', 2, 'Location', 'northeast', 'FontSize', 7);
lgd.Title.String = 'Hs (Sz width)';
%lgd.FontSize = textsz - 4;
text(-0.8, 0.375, {'cBathy', '(hErr < 0.5)'}, 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', textsz)
set(gca, 'FontSize', textsz)
ylim([0 0.75])
xlim([-0.25 1.5])
xticks([0:0.5:1])
yticks([0 0.5])
grid on
xticklabels({''})
text(-0.22, 0.65, 'a', 'FontSize', textsz+5)
    
ax = subplot(423); ax.Position = [0.17; 0.56; 0.4; 0.2]
h=bar(x,hnlin,'stacked'); for ii = 1:6;h(ii).FaceColor = c(ii,:);end
text(-0.8, 0.375, {'cBathy', '(nonlinear)'}, 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', textsz)
set(gca, 'FontSize', textsz)
ylim([0 0.75])
xlim([-0.25 1.5])
xticks([0:0.5:1])
yticks([0 0.5])
grid on
xticklabels({''})
text(-0.22, 0.65, 'c', 'FontSize', textsz+5)

ax = subplot(425); ax.Position = [0.17; 0.35; 0.4; 0.2]
h=bar(x,hbp,'stacked'); for ii = 1:6;h(ii).FaceColor = c(ii,:);end
text(-0.8, 0.375, 'cBathyCT', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', textsz)
text(0.625, -0.23, 'RMSE (m)', 'HorizontalAlignment', 'center', 'FontSize', textsz)
set(gca, 'FontSize', textsz)
ylim([0 0.75])
xlim([-0.25 1.5])
xticks([0:0.5:1])
yticks([0 0.5])
grid on
xtickangle(0)
text(-0.22, 0.65, 'e', 'FontSize', textsz+5)

ax = subplot(427); ax.Position = [0.17; 0.06; 0.4; 0.2]
h=bar(xb,bias,'stacked'); for ii = 1:6;h(ii).FaceColor = c(ii,:);end
text(-1.3, 0.375, 'cBathyCT', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', textsz)
text(0.125, -0.24, 'Bias (m)', 'HorizontalAlignment', 'center', 'FontSize', textsz)
set(gca, 'FontSize', textsz)
ylim([0 0.75])
xlim([-0.75 1])
xticks([-0.5:0.5:1])
yticks([0 0.5])
grid on
xtickangle(0)
text(-0.72, 0.65, 'g', 'FontSize', textsz+5)

    
ax = subplot(422); ax.Position =  [0.59; 0.77; 0.4; 0.2]
h=bar(x,hcb3,'stacked'); for ii = 1:6;h(ii).FaceColor = c(ii,:);end
tl = title('3 Flights', 'FontSize', textsz);
set(gca, 'FontSize', textsz)
ylim([0 0.75])
xlim([-0.25 1.5])
xticks([0:0.5:1])
yticks([0 0.5])
grid on
xticklabels({''})
yticklabels({''})
text(-0.22, 0.65, 'b', 'FontSize', textsz+5)
    
ax = subplot(424); ax.Position = [0.59; 0.56; 0.4; 0.2]
h=bar(x,hnlin3,'stacked'); for ii = 1:6;h(ii).FaceColor = c(ii,:);end
set(gca, 'FontSize', textsz)
ylim([0 0.75])
xlim([-0.25 1.5])
xticks([0:0.5:1])
yticks([0 0.5])
grid on
xticklabels({''})
yticklabels({''})
text(-0.22, 0.65, 'd', 'FontSize', textsz+5)

ax = subplot(426); ax.Position = [0.59; 0.35; 0.4; 0.2]
h=bar(x,hbp3,'stacked'); for ii = 1:6;h(ii).FaceColor = c(ii,:);end
text(0.625, -0.23, 'RMSE (m)', 'HorizontalAlignment', 'center', 'FontSize', textsz)
set(gca, 'FontSize', textsz)
ylim([0 0.75])
xlim([-0.25 1.5])
xticks([0:0.5:1])
yticks([0 0.5])
grid on
yticklabels({''})
xtickangle(0)
text(-0.22, 0.65, 'f', 'FontSize', textsz+5)

ax = subplot(428); ax.Position = [0.59; 0.06; 0.4; 0.2]
h=bar(xb,bias3,'stacked'); for ii = 1:6;h(ii).FaceColor = c(ii,:);end
text(0.125, -0.27, 'Bias (m)', 'HorizontalAlignment', 'center', 'FontSize', textsz)
set(gca, 'FontSize', textsz)
ylim([0 0.75])
xlim([-0.75 1])
xticks([-0.5:0.5:1])
yticks([0 0.5])
grid on
yticklabels({''})
xtickangle(0)
text(-0.72, 0.65, 'h', 'FontSize', textsz+5)


set(0,'DefaultFigureColor','remove')
print('/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/Figures/FIG8.pdf','-dpdf', '-r300', '-bestfit')
    
%close all

%% Fig 10 - Gamma Sensitivity

[hFig, ax1] = makeFig(7,5,3,2)

ax = subplot(3,2,[1 3 5]); ax.Position = [0.1; 0.05; 0.38; 0.87]
load('/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/Video_gamma.mat')
  
rgbp = sort(RMSE_gamma_bp(:,id_Torrey_gamma)');
n=length(rgbp);
    gamma = gamma_opt;
    %patch([gamma fliplr(gamma)]', [nanstd(RMSE_gamma_bp(:,id_Torrey_gamma),[],2)+nanmean(RMSE_gamma_bp(:,id_Torrey_gamma),2); flipud(nanmean(RMSE_gamma_bp(:,id_Torrey_gamma),2)-nanstd(RMSE_gamma_bp(:,id_Torrey_gamma),[],2))],[0.6 0.6 0.6], 'EdgeColor', 'none');alpha 0.4
    hold on
    %plot(gamma, nanmedian(RMSE_gamma_bp(:,id_Torrey_gamma),2), 'b.-', 'LineWidth', 3)
%     plot(gamma, rgbp(round(n*.5 - n*.25),:), 'k.-', 'LineWidth', 3)
%     plot(gamma, rgbp(round(n*.5 + n*.25),:), 'k.-', 'LineWidth', 3)
%     plot(gamma, rgbp(round(n*.5 - n*.45),:), 'r.-', 'LineWidth', 3)
%     plot(gamma, rgbp(round(n*.5 + n*.45),:), 'r.-', 'LineWidth', 3)
% patch([gamma fliplr(gamma)]', [rgbp(round(n*.5 + n*.45),:) fliplr(rgbp(round(n*.5 - n*.45),:))], colors(3,:), 'EdgeColor', colors(3,:));alpha 0.4
% patch([gamma fliplr(gamma)]', [rgbp(round(n*.5 + n*.25),:) fliplr(rgbp(round(n*.5 - n*.25),:))], [0.8 0.8 0.8], 'EdgeColor', [0.8 0.8 0.8]);alpha 0.8
% plot(gamma, rgbp(round(n*.5),:), '.-', 'Color', colors(5,:), 'LineWidth', 3)
 
%     for ii = 1:length(gamma)
%         patch([gamma(ii)-0.01 gamma(ii)+0.01 gamma(ii)+0.01 gamma(ii)-0.01]', [rgbp(round(n*.5 + n*.45),ii) rgbp(round(n*.5 + n*.45),ii) rgbp(round(n*.5 + n*.25),ii) rgbp(round(n*.5 + n*.25),ii)], colors(3,:), 'EdgeColor', colors(3,:));alpha 0.6
%         patch([gamma(ii)-0.01 gamma(ii)+0.01 gamma(ii)+0.01 gamma(ii)-0.01]', [rgbp(round(n*.5 - n*.45),ii) rgbp(round(n*.5 - n*.45),ii) rgbp(round(n*.5 - n*.25),ii) rgbp(round(n*.5 - n*.25),ii)], colors(3,:), 'EdgeColor', colors(3,:));alpha 0.6
%         patch([gamma(ii)-0.01 gamma(ii)+0.01 gamma(ii)+0.01 gamma(ii)-0.01]', [rgbp(round(n*.5 + n*.25),ii) rgbp(round(n*.5 + n*.25),ii) rgbp(round(n*.5 - n*.25),ii) rgbp(round(n*.5 - n*.25),ii)], [0.5 0.5 0.5], 'EdgeColor', [0.8 0.8 0.8]);alpha 0.6
%         plot([gamma(ii)-0.01 gamma(ii)+0.01], [rgbp(round(n*.5),ii) rgbp(round(n*.5),ii)], 'Color', colors(5,:), 'LineWidth', 3)
%     end

    bx = boxplot(RMSE_gamma_bp([1:2:9 12:2:22], id_Torrey_gamma)', gamma([1:2:9 12:2:22]),'OutlierSize',0.00000000001)
    hold on
    set(findobj(bx,'Tag','Box'),'LineWidth', 1.5)
    set(findobj(bx,'Tag','Lower Whisker'),'LineWidth', 1.5, 'LineStyle', '-')
    set(findobj(bx,'Tag','Upper Whisker'),'LineWidth', 1.5, 'LineStyle', '-')
    plot(nanmedian(RMSE_gamma_bp([1:2:9 12:2:22],id_Torrey_gamma),2), '.-', 'Color', colors(5,:), 'LineWidth', 4)
    %xtickslabels(num2str([gamma([1:2:9 12:2:22])]))
    xlabel('\gamma')
    yl = ylabel('$\langle$ RMSE $\rangle$ (m)', 'Interpreter', 'latex');
    p=get(yl, 'Position');p(1) = p(1)-1; p(2) = 0.5;
    set(yl, 'Position', p)

    set(gca, 'FontSize', textsz)
    ylim([0 1])

    tl = title('RMSE Sensitivity to \gamma in BP', 'FontSize', textsz)
    tl.Position(2) = 1.035
    grid on
    text(10.5, 0.96, 'a', 'FontSize', textsz + 5)

    load('/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/DATA/lidar_gamma_new.mat')
    clear id_* date p
    for rr = 1:length(hover)
        date(rr) = str2double(hover(rr).name(1:8));

    end
    id_20210712 = find(date == 20210712);
    id_20211026 = find(date == 20211026);
    id_20211102 = find(date == 20211102);
    id_20211116 = find(date == 20211116);
    id_20211202 = find(date == 20211202);
    id_20211215 = find(date == 20211215);


    c = turbo(7);c(1,:)=[];c(end,:)=[];
    ax = subplot(322); ax.Position = [0.61; 0.675; 0.37; 0.25]
    hold on
    id_run = id_20211102;
    x=repmat([-155:0.1:0],length(id_run), 1);
    gamma=NaN(size(x));
    for ii = id_run
        if hover(ii).SWL > 0.25
        for jj = 1:size(x,2)-id_run(1)
            if any(x(ii-id_run(1)+1,jj)== hover(ii).x)
                id=find(x(ii-id_run(1)+1,jj)== hover(ii).x);
                gamma(ii-id_run(1)+1,jj) = hover(ii).gamma(id);
            end
        end
        end
    end
    
    aa=[hover(id_run).SWL]; if aa(1) > aa(end);id_run = fliplr(id_run);end
    jj=0;
    for ii = id_run
        if hover(ii).SWL>0.25
            jj=jj+1;
            p(jj) = scatter(-hover(ii).x, hover(ii).gamma, 7, c(jj,:),'filled');%, hover(ii).SWL*ones(length(hover(ii).h),1))
        end
    end
    
    mean_transect = nanmean(gamma,1); id_cutoff = 951;
    gamma_mean = nanmean(mean_transect(1:id_cutoff),2);
    min(min(gamma(:, 1:id_cutoff)))
    max(max(gamma(:, 1:id_cutoff)))
    plot([60 125], [gamma_mean gamma_mean], 'k', 'LineWidth', 2)
    text(75, 0.41,texlabel(sprintf('gamma = %.2f', gamma_mean)),'FontSize', textsz)
    %title('Hs = 1.12m')
    %lgd = legend(p,string(round([hover(id_run).SWL],2)), 'NumColumns', 2, 'box', 'off');
    %lgd.Title.String = 'SWL';
    xticks([40:20:120])
    xticklabels({''})
    ylim([0.2 1])
    xlim([30 140])
    set(gca, 'FontSize', textsz)
    yl = ylabel('\gamma', 'FontSize', textsz);
    p=get(yl, 'Position');p(1) = p(1)-10;
    set(yl, 'Position', p)
    grid on
    
    tl = title('LiDAR observed \gamma', 'FontSize', textsz)
    tl.Position(2) = 1.09
    text(130, 0.9, 'b', 'FontSize', textsz + 5)
    set(ax, 'box', 'on')
    
    
    ax = subplot(324);  ax.Position = [0.61; 0.4045; 0.37; 0.25]
    hold on
    id_run = id_20211215;
    x=repmat([-155:0.1:0],length(id_run), 1);
    gamma=NaN(size(x));
    for ii = id_run
        if hover(ii).SWL > 0.25
        for jj = 1:size(x,2)-id_run(1)
            if any(x(ii-id_run(1)+1,jj)== hover(ii).x)
                id=find(x(ii-id_run(1)+1,jj)== hover(ii).x);
                gamma(ii-id_run(1)+1,jj) = hover(ii).gamma(id);
            end
        end
        end
    end
    
    aa=[hover(id_run).SWL]; if aa(1) > aa(end);id_run = fliplr(id_run);end
    jj=0;
    for ii = id_run
        if hover(ii).SWL>0.25
            jj=jj+1;
            p(jj) = scatter(-hover(ii).x, hover(ii).gamma, 7, c(jj,:), 'filled');%, hover(ii).SWL*ones(length(hover(ii).h),1))
        end
    end
    
    mean_transect = nanmean(gamma,1); id_cutoff = 951;
    gamma_mean = nanmean(mean_transect(1:id_cutoff),2);
    min(min(gamma(:, 1:id_cutoff)))
    max(max(gamma(:, 1:id_cutoff)))
    plot([60 140], [gamma_mean gamma_mean], 'k', 'LineWidth', 2)
    text(62, 0.32,texlabel(sprintf('gamma = %.2f', gamma_mean)),'FontSize', textsz)
    %title('Hs = 1.92m')
    %lgd = legend(p,string(round([hover(id_run).SWL],2)), 'NumColumns', 2, 'box', 'off');
    %lgd.Title.String = 'SWL';
    xticks([40:20:120])
    xticklabels({''})
    ylim([0.2 1])
    xlim([30 140])
    set(gca, 'FontSize', textsz)
    yl = ylabel('\gamma', 'FontSize', textsz);
    p=get(yl, 'Position');p(1) = p(1)-10;
    set(yl, 'Position', p)
    grid on
    text(130, 0.9, 'c', 'FontSize', textsz + 5)
    set(ax, 'box', 'on')
    

    ax = subplot(326);  ax.Position = [0.61; 0.134; 0.37; 0.25]
    hold on
    id_run = id_20211026;
    x=repmat([-155:0.1:0],length(id_run), 1);
    gamma=NaN(size(x));
    for ii = id_run
        if hover(ii).SWL > 0.25
        for jj = 1:size(x,2)-id_run(1)
            if any(x(ii-id_run(1)+1,jj)== hover(ii).x)
                id=find(x(ii-id_run(1)+1,jj)== hover(ii).x);
                gamma(ii-id_run(1)+1,jj) = hover(ii).gamma(id);
            end
        end
        end
    end
    
    aa=[hover(id_run).SWL]; if aa(1) > aa(end);id_run = fliplr(id_run);end
    jj=0;
    for ii = id_run
        if hover(ii).SWL>0.25
            jj=jj+1;
            p(jj) = scatter(-hover(ii).x, hover(ii).gamma, 7, c(jj,:), 'filled');%, hover(ii).SWL*ones(length(hover(ii).h),1))
        end
    end
    
    mean_transect = nanmean(gamma,1); id_cutoff = 951;
    gamma_mean = nanmean(mean_transect(1:id_cutoff),2);
    min(min(gamma(:, 1:id_cutoff)))
    max(max(gamma(:, 1:id_cutoff)))
    plot([60 140], [gamma_mean gamma_mean], 'k', 'LineWidth', 2)
    text(62, 0.32,texlabel(sprintf('gamma = %.2f', gamma_mean)),'FontSize', textsz)
    %title('Hs = 2.15m')
    %lgd = legend(p,string(round([hover(id_run).SWL],2)), 'NumColumns', 2, 'box', 'off');
    %lgd.Title.String = 'SWL';
    xlabel('Cross-shore Distance (m)')
    ylim([0.2 1])
    xlim([30 140])
    set(gca, 'FontSize', textsz)
    yl = ylabel('\gamma', 'FontSize', textsz);
    p=get(yl, 'Position');p(1) = p(1)-10;
    set(yl, 'Position', p)
    grid on
    text(130, 0.9, 'd', 'FontSize', textsz + 5)
    xticks([40:20:120])
    set(ax, 'box', 'on')
    

    set(0,'DefaultFigureColor','remove')
    print('/Volumes/LANGE_Passport/Remote_surfzone_bathy_paper/Figures/FIG10.pdf','-dpdf', '-r300', '-bestfit')
    
    close all
%% Appendix Breakpoint Method
%breakpoint_method_model_plot


%% Functions
function [rmse, skill, bias] = error(obs, pred)
    rmse = sqrt(sum((pred - obs).^2)/length(obs));
    skill = 1 - sum((pred - obs).^2)./sum((obs - mean(obs)).^2);
    bias = sum(pred - obs) / length(obs);
end
function [r,s,b] = calc_errors(obs, pred, lim)
    tempx = pred;
    tempx = tempx(lim(1):lim(2));
    ax = obs;
    ax = ax(lim(1):lim(2));
    ax(find(isnan(tempx)))=[];
    tempx(find(isnan(tempx)))=[];
    tempx(find(isnan(ax)))=[];
    ax(find(isnan(ax)))=[];
    [r,s,b]=error(ax, tempx);
end
