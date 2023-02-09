%% Compute bathymetry inversion from cBathy and timestack wave crests
%% Input
%   survey (structure)
%       date - YYYYMMDD
%       moplines - [579, 580, 581, 582] - should be the same naming
%       convention as timestack images and cBathy grid x
%       z - [moplines, timestack cross-shore] - should be interpolated onto
%       timestack cross-shore grid (0.1m spacing)
%
%   cBathy (structure)
%       date - YYYYMMDD
%       hover - hover number of the day (used to match with timestack name)
%       x - x coordinates from localX [alongshore coordinates, cross-shore coordinates](e.g., [775:5:-85])
%       y - y coordinates from localY (MOP numbers) [alongshore coordinates, cross-shore coordinates] (e.g., [576:0.2:588.8])
%       z - depth values (-bathy.fCombined.h), with tidal correction [alongshore coordinates, cross-shore coordinates]
%       zerr - 95per confidence interval (bathy.fCombined.hErr) [alongshore coordinates, cross-shore coordinates]
%       tide - tide level during hover (pulled from cBathy stucture array)
%
%   timestacks - located in [date_dir '/timestacks/data]
%
%
%% Installation
%
%   setup_nctoolbox
%   NOAA tide gauge number (line 154)
%   search box size for beginning of wave crests - higher frequency days may require smaller size (line 216)
%   minimum time duration for wave crest (line 287)
%   5m Gaussian smoothing (line 431)
%
%
%% Output
%   Video_bathy (structure)
%       date - YYYYMMDD
%       location - ie. 'Torrey' (pulled from timestack name)
%       flight - flight number
%       mop - mop number
%       x10 - [0:0.1:500]m
%       survey - z - Depth pulled from survey on MOP line
%       cbathy - 
%       tide - tide level (pulled from cbathy)
%       min_tide - minimum tidal level of the day (for subaerial survey offshore location)
%       crests - 
%           t - time for wave tracks (sec)
%           x - cross-shore location for wave tracks (m) (0 is offshore, 500 is onshore)
%           c_15m_avg - phase speed of wave tracks - interp to x10 (from RANSAC slope over 15m distance)
%       bp - breakpoint index of wave tracks
%       xshift - cross-shore shift index to match subaerial survey with subaqueous bathymetry
%       h_avg - interped to x10
%           lin - linear crest-tracking c = sqrt(gh)
%           nlin - nonlinear crest-tracking c = sqrt(gh(1+0.42))
%           bp - breakpoint transition crest-tracking c = sqrt(gh(1+gamma(x)))
%       gamma - transition between 0 and 0.42 for wave tracks (interped to x10)
%       gamma_mean - gamma(x) - mean of step function gamma for wave tracks
%       composite - constructed bathymetry with subaerial survey with ...
%           cbathy_hErr - cBathy with hErr < 0.5m 
%           cbathy_gamma - cBathy with breaking region removed
%           cbathy_nlin - cBathy_gamma with gamma(x) correction
%           cbathyCT - breakpoint transition crest-tracking surfzone bathymetry and cBathy offshore
%       lims - index of [1st BP valid onshore point, onshore cutoff of breaking, offshore cutoff of breaking]
%       Error - 
%           RMSE - root-mean-square error
%           Skill - skill of estimated versus observed survey bathymetry
%           Bias - bias of estimated bathymetry
%
%           insz - inner surfzone (between shoreline and end of active wave breaking - gamma(x)=0.42)
%           break - breaking region (active wave breaking region - where gamma changes)
%           sz - full surfzone region (between shoreline and wave breaking region)
%           full - [0 500] region
%           offshore - offshore region (between beginning of wave breaking and 500m - gamma(x) = 0)
%           
%           cb - default cBathy (no region removed)
%           cb_hErr - cBathy with hErr > 0.5m removed     
%           cb_gamma - cBathy with breaking region removed
%           lin - linear crest-tracking
%           nlin - nonlinear crest-tracking
%           bp - breakpoint transition crest-tracking
%           comp_hErr - composite cBathy_hErr
%           comp_gamma - composite cBathy_gamma
%           comp_nlin - composite cBathy_nlin
%           comp_CT - composite BP + cBathy
%
%% Copyright
% (c) Lange, Athina M.Z., Fiedler, Julia W., Merrifield, Mark A., Guza, R.T. (2022)
% UAV video-based estimates of nearshore bathymetry, submitted to Coastal Engineering
%% Determine number of timestacks/bathys to process
% find all timestacks that have been processed
close all
clearvars -except stationStr
% define local path
local_dir = uigetdir(pwd,'Local Directory');
cd(local_dir)

% define location of data and codes
code_dir = fullfile(local_dir, 'CODES');
addpath(genpath(code_dir))
data_dir = fullfile(local_dir, 'DATA');
addpath(genpath(data_dir))


files = dir(fullfile(data_dir, 'timestacks', 'processed'));

% remove empty files
for nn = length(files):-1:1 
    if files(nn).name(1) == '.'
        files(nn)=[];
    end
end
Video = struct('date',repmat({''},length(files),1));

%setup_nctoolbox
%% get pixel resolution
dx=0.5;
[cutoff]  = pixel_res(files, data_dir, local_dir, dx);
%% Wave crest detection
for rr = 1:length(files)
    tic
    rr
    % save time and overwriting
    close all
    clearvars -except *_dir Video files rr flight_correction survey cbathy min_tide cutoff tide
    
    % determine date / location / flight number / mop line
    ab = split(files(rr).name, '_');
    date = str2num(ab{1});
    location = ab{2};
    flight = str2num(ab{3});
    mop = str2num(ab{4});
    
    % save to Video file
    Video(rr).date = date;
    Video(rr).location = location;
    Video(rr).flight = flight;
    Video(rr).mop = mop;

    %% Bathy from survey
    Video(rr).x10 = [0:0.1:500]';
    Video(rr).survey.z = survey(find(date == [survey.date])).z(survey(find(date == [survey.date])).moplines==mop,:)';
    
    %% cBathy 
    id=find([cbathy.date]==Video(rr).date & [cbathy.hover]==Video(rr).flight);
    Video(rr).cbathy.z = interp1(cbathy(id).x(find(round([cbathy(id).y(:,1)],1)==round(Video(rr).mop,1)),:),cbathy(id).z(find(round([cbathy(id).y(:,1)],1)==round(Video(rr).mop,1)),:), Video(rr).x10);
    Video(rr).cbathy.zerr = interp1(cbathy(id).x(find(round([cbathy(id).y(:,1)],1)==round(Video(rr).mop,1)),:),cbathy(id).zerr(find(round([cbathy(id).y(:,1)],1)==round(Video(rr).mop,1)),:), Video(rr).x10);
    Video(rr).tide = cbathy(id).tide;
     
    % find lowest tide of day ~ when survey was taken to provide
    % subaerial survey limit
    tide_gauge = '9410230'; % NOAA tide gauge at Scripps Institution of Oceanography
    aa=char(string(date));
    [~,~,verified,~,~] = getNOAAtide(datetime(str2num(aa(1:4)), str2num(aa(5:6)), str2num(aa(7:8)), 8,0,0), datetime(str2num(aa(1:4)), str2num(aa(5:6)), str2num(aa(7:8)), 18,0,0),tide_gauge);
    min_tide = min(verified);
    clear aa
    Video(rr).min_tide=min_tide;
        
    %% Compute wave celerity from timestack wave crests
    %%% Load timestack
    if Video(rr).flight < 10
        gt = imread(fullfile(data_dir, 'timestacks','processed', [char(string(Video(rr).date)) '_' char(Video(rr).location) '_0' char(string(Video(rr).flight)) '_' char(string(Video(rr).mop)) '_prediction.jpg']));
        [fi,fcmap] = imread(fullfile(data_dir, 'timestacks','data', [char(string(Video(rr).date)) '_' char(Video(rr).location) '_0' char(string(Video(rr).flight)) '_' char(string(Video(rr).mop)) '.png']));
    else
        gt = imread(fullfile(data_dir, 'timestacks','processed', [char(string(Video(rr).date)) '_' char(Video(rr).location) '_' char(string(Video(rr).flight)) '_' char(string(Video(rr).mop)) '_prediction.jpg']));
        [fi,fcmap] = imread(fullfile(data_dir, 'timestacks','data', [char(string(Video(rr).date)) '_' char(Video(rr).location) '_' char(string(Video(rr).flight)) '_' char(string(Video(rr).mop)) '.png']));
    end   
    % wave crests should be a binary image
    if size(gt,3) ~= 1
        gt = rgb2gray(gt); 
    end
    % check that all png 24bit
    if ~isempty(fcmap)
        fi = ind2rgb(fi,fcmap);
    end
                
    Video(rr).timestack = fi;
    Video(rr).waves = gt;
    
    %%% Find the average line of the wave crests as determined by the neural network
    gt = Video(rr).waves;
    % binarize and skeletonize the wave crests
    gt_clean = bwskel(imbinarize(imcomplement(gt)),'MinBranchLength',20);
    
    % turn the wave crests into 0 (wave crest) or 255 (not)
    gt_gray = ones(size(gt))*255;
    gt_gray(gt_clean==1)=0;
    
    %%% Find midpoint of each x segment - only want 1 value for every x point
    for tt = 1:size(gt,2)
        for xx = 1:size(gt,1)-20
            if gt_gray(xx,tt)==0 % if marked as wave crest
                aa = find(gt_gray(xx:xx+20,tt)==0);
                if length(aa)==1 
                    midpt = 0;
                    aa=[];
                elseif length(aa)==2
                    aa(1)=[];
                elseif mod(length(aa),2)==1 % odd
                    midpt = aa((length(aa)-1)/2);
                    aa(midpt+1==aa)=[];
                elseif mod(length(aa),2)==0 % even
                    midpt = aa((length(aa))/2)-1;
                    aa(midpt==aa)=[];
                end
                 gt_gray(xx+aa-1,tt)=255; % find which points should be unmarked
            end
        end
    end


    %%% Find peaks - 1st point of wave crest
    clear peaks_t peaks_x
    box = 20; % can choose to change search box size HERE - higher frequency days may require smaller size 
    
    % go through full image and stop when you find a 0 pixel (wave crest)
    % if all the pixels previously in time and space within the search
    % box are 255, then it is the beginning of a crest, and save the
    % coordinates. 
    for tt = box+1:size(gt,2)
        for xx = box+1:size(gt,1)
            if gt_gray(xx,tt)==0
                if (sum(gt_gray(xx-box:xx, tt-box:tt), 'all')+255)/((box+1)^2) == 255 
                    % any other 0 value will skew results to not be = 255
                    peaks_t(tt) = tt;
                    peaks_x(tt) = xx;
                    tt = tt+1;
                end
            end
        end
    end
    % remove any weirdness
    peaks_x(peaks_t==0)=[];
    peaks_t(peaks_t==0)=[];


    %%% Find wave trains from peaks - change box size for next wave point HERE
    clear crests_t crests_x
    
    for nn = 1:length(peaks_t) % loop through all peaks found
        
        ll=1;
        crests_t(ll,nn) = peaks_t(nn);
        crests_x(ll,nn) = peaks_x(nn);
    
        D=1;
        while min(D) < 20 % loop through wave train
            tt = crests_t(ll,nn);
            xx = crests_x(ll,nn);
            ll=ll+1;
            % if the wave train is at the right most edge of the image,
            % stop looking - avoids throwing an error here
            if xx > size(gt_gray,1)-30 | tt > size(gt_gray,2)-5 
                break
            end
            % within a box size onshore and in future time, find all
            % wave crest pixels. 
            [idx,idt]=find(gt_gray(xx+1:xx+30, tt+1:tt+5)==0);
            if length(idx) < 2 
                D=51;
                break
            end
            X = [1,1];
            Y = [idt,idx];%Y(1,:)=[];
            % find the point that has the shortest euclidean distance
            % from the previously found point 
            [D,~] = pdist2(X,Y,'euclidean','Smallest',1);
            if length(find(min(D)==D))>1
                id = find(min(D)==D); id=id(1);
                crests_t(ll,nn) = tt+idt(id);
                crests_x(ll,nn) = xx+idx(id);
            else
                crests_t(ll,nn) = tt+idt(find(min(D)==D));
                crests_x(ll,nn) = xx+idx(find(min(D)==D));
            end
        end
    end
    % remove any weirdness
    crests_t(crests_t==0)=NaN;
    crests_x(crests_x==0)=NaN;


    %%% Clean up short segements - change length of segment HERE
    for nn = length(peaks_t):-1:1
        if find(~isnan(crests_t(:,nn))) < 75
            crests_t(:,nn)=[];
            crests_x(:,nn)=[];
        end
    end


    %%% Get velocity
    
    % need to go from pixel units to real world units (0.1m, 0.1s)
    crests_x = crests_x./10;
    crests_t = crests_t./10;
    clear c_new c_short
    % compute velocity using RANSAC
    c_new = NaN(size(crests_t,1)-1, size(crests_t,2));
    for nn = 1:size(crests_t,2)
        for ll = 1:length(crests_x(:,nn))-1
            if isnan(crests_x(ll,nn)) == 1
                break
            else
                clear id
                id = find(crests_x(:,nn) >= crests_x(ll,nn)-7.5 & crests_x(:,nn) <= crests_x(ll,nn)+7.5);
                if ll == 1 
                    c_new(ll,nn) = diff(crests_x(ll:ll+1,nn))./diff(crests_t(ll:ll+1,nn));
                elseif isnan(crests_x(ll+1,nn)) == 1
                    break
                elseif length(id) > 1
                    points = [crests_x(id,nn) crests_t(id,nn)];
                    sampleSize = 2; % number of points to sample per trial
                    maxDistance = 0.1; % max allowable distance for inliers
    
                    fitLineFcn = @(points) polyfit(points(:,1),points(:,2),1); % fit function using polyfit
                    evalLineFcn = @(model, points) sum((points(:, 2) - polyval(model, points(:,1))).^2,2);
    
                    [modelRANSAC, inlierIdx] = ransac(points,fitLineFcn,evalLineFcn, sampleSize,maxDistance);
                    
                    c_new(ll,nn)= 1/modelRANSAC(1);
                
                else
                    c_new(ll,nn)=diff(crests_x(ll:ll+1,nn))./diff(crests_t(ll:ll+1,nn));
                end
            end
        end
    end
    for nn = 1:size(crests_t,2)
        x10 = Video(rr).x10;
        aa = interp1(crests_x(~isnan(c_new(:,nn)),nn), c_new(~isnan(c_new(:,nn)),nn), x10, 'nearest');
        c_short(:,nn) = fliplr(aa);
    end
    
    % images are named from top left corner, but we define x10 going from 
    % bottom to top (shore -> offshore)
    x10 = (flipud(x10))'; 

    %%% find cutoffs where not enough data present - change % for cutoff HERE
    for nn = 1:length(x10); counts(nn) = sum(~isnan(c_short(nn,:))); end
    max_id = round(median(find(max(counts)==counts)));
    cutoff = max(counts).*.10; ab=abs(counts-max(counts).*.10); 
    min_diff = min(ab(max_id:end)); for ii = max_id:1:length(x10); if ab(ii) == min_diff; offshore_cutoff = ii; break; end; end
    min_diff = min(ab(1:max_id)); for ii = max_id:-1:1; if ab(ii) == min_diff; onshore_cutoff = ii; break; end; end
    c_short(1:onshore_cutoff,:)= NaN;
    c_short(offshore_cutoff:end,:)=NaN;  

    Video(rr).crests.t = crests_t;
    Video(rr).crests.x = crests_x;
    Video(rr).crests.t_for_c = t_slots;
    Video(rr).crests.x_for_c = repmat(x10', 1,size(crests_t,2));
    Video(rr).crests.c_15m_avg = c_short;


    %% Breakpoint Finder
    bp = breakpt_calculator(Video, rr);
    Video(rr).bp = bp;
    % Blue dot: Unbroken, dot at end of track
    % Green dot: Breakpoint location along wave track
    % Red dot: Entire wave track broken, dot at beginning of track
    
    toc
end

clear Video_short
Video_short = Video;
Video_short = rmfield(Video_short, 'timestack');
Video_short = rmfield(Video_short, 'waves');
clear Video
Video_base = Video_short;

save(fullfile(data_dir, 'Video_base.mat'), 'Video_base', '-v7.3')

%% Surfzone inversion and composite bathy
load(fullfile(data_dir, 'cutoff_pixres.mat'))
load(fullfile(data_dir, 'Video_base.mat'))
Video = Video_base;
for rr = 1:length(Video)
    %%
    rr
    tic
    % save time and overwriting
    close all
    clearvars -except *_dir Video* files rr flight_correction survey cbathy min_tide cutoff
    
    date = Video(rr).date;
    mop = Video(rr).mop;
    location = Video(rr).location;
    %%% Cutoff offshore based on pixel resolution
    id = find([cutoff.date] == Video(rr).date & [cutoff.hover] == Video(rr).flight);
    id_mop = find(cutoff(id).mopgrid == round(Video(rr).mop));
       
    min_tide=Video(rr).min_tide;
       
    gamma_H = 0.42; % gamma = H/h value 
    x10 = Video(rr).x10;
    % images are named from top left corner, but we define x10 going from 
    % bottom to top (shore -> offshore)
    x10 = (flipud(x10))';
    %% Surfzone depth inversion - default c
    %%% compute depth inversion with $c = \sqrt(gh(1+\gamma))$
    c_real = Video(rr).crests.c_15m_avg;
    break_loc = (Video(rr).bp);
    % set completely broken to 0
    if ~isempty(find(isnan(break_loc)==1))
        break_loc(find(isnan(break_loc)==1)) = 0;
    end
    
    for nn = 1:size(c_real,2)
        % compute array of 0 or 1.42 depending on if broken or not
        corr = (1+gamma_H)*ones(size(x10))';
        % if not all broken
        if break_loc(nn) ~=0
            corr(1:break_loc(nn))=1;
        end

        % Get bathy inversion
        h_crest_lin(:,nn) = c_real(:,nn).^2./(9.81);
        h_crest_lin(find(h_crest_lin(:,nn)>10),nn)=NaN;

        h_crest_nlin(:,nn) = c_real(:,nn).^2./(9.81.*(1+gamma_H));
        h_crest_nlin(find(h_crest_nlin(:,nn)>10),nn)=NaN;
     
        h_crest_bp(:,nn) = c_real(:,nn).^2./(9.81.*corr);
        h_crest_bp(find(h_crest_bp(:,nn)>10),nn)=NaN;
    
    end
        
    %%% Compute h avg 
    af=find(~isnan(nanmean(h_crest_lin,2)));
    smooth_factor = 5*10;
    % 5m Gaussian smoothing after averaging
    h_crest_lin_avg = smoothdata(nanmean(h_crest_lin,2), 'gaussian', smooth_factor);
    h_crest_nlin_avg = smoothdata(nanmean(h_crest_nlin,2), 'gaussian', smooth_factor);
    h_crest_bp_avg = smoothdata(nanmean(h_crest_bp,2), 'gaussian', smooth_factor);
    
    % remove artefacts of smoothing
    h_crest_lin_avg(1:af(1))=NaN;
    h_crest_nlin_avg(1:af(1))=NaN;
    h_crest_bp_avg(1:af(1))=NaN;
    h_crest_lin_avg(af(end):end)=NaN;
    h_crest_nlin_avg(af(end):end)=NaN;
    h_crest_bp_avg(af(end):end)=NaN;
    
    % flip to match x10 in structure
    if isrow(h_crest_lin_avg); h_crest_lin_avg = h_crest_lin_avg'; end
    if isrow(h_crest_nlin_avg); h_crest_nlin_avg = h_crest_nlin_avg'; end
    if isrow(h_crest_bp_avg); h_crest_bp_avg = h_crest_bp_avg'; end

    % go from depth to NAVD88m with tidal correction
    h_crest_lin_avg = flipud(-h_crest_lin_avg) + Video(rr).tide;
    h_crest_nlin_avg = flipud(-h_crest_nlin_avg) + Video(rr).tide;
    h_crest_bp_avg = flipud(-h_crest_bp_avg) + Video(rr).tide;

    % remove any downsloping curves
    [~,i]=max(h_crest_lin_avg); h_crest_lin_avg(1:i-1)=NaN;
    [~,i]=max(h_crest_nlin_avg); h_crest_nlin_avg(1:i-1)=NaN;
    [~,i]=max(h_crest_bp_avg); h_crest_bp_avg(1:i-1)=NaN;

    if ~isnan(cutoff(id).cbathy(id_mop))
        Video(rr).cbathy.z(find(Video(rr).x10 == -cutoff(id).cbathy(id_mop)):end) = NaN;
    else
        Video(rr).cbathy.z(1:end) = NaN;
    end
  
    %%% x-shift to match subaerial surveys
    id = find(~isnan(h_crest_bp_avg)==1); % find first non-nan of predicted bathy
    if isempty(id)
        xshift = NaN;
    elseif h_crest_bp_avg(id(1)) < min_tide % if predicted bathy doesn't reach to min tide
        [~,id2]=min(abs((h_crest_bp_avg(id(1)))-Video(rr).survey.z));
        xshift =-(id2-id(1));
    else % when predicted bathy does reach min tide
        [~,ab]=min(abs(h_crest_bp_avg-min_tide));
        [~,ac]=min(abs(Video(rr).survey.z-min_tide));
        xshift = ab-ac;
        clear ab ac
    end
    Video(rr).xshift = xshift;
    
    if xshift < 0 % get rid of shoreward points
        xshift = abs(xshift);
        Video(rr).h_avg.lin = [NaN(xshift,1); h_crest_lin_avg(1:end-xshift)];
        Video(rr).h_avg.nlin = [NaN(xshift,1); h_crest_nlin_avg(1:end-xshift)];
        Video(rr).h_avg.bp = [NaN(xshift,1); h_crest_bp_avg(1:end-xshift)];
    elseif xshift == 0
        Video(rr).h_avg.lin = h_crest_lin_avg;
        Video(rr).h_avg.nlin = h_crest_nlin_avg;
        Video(rr).h_avg.bp = h_crest_bp_avg;
    elseif xshift > 0 % add shoreward points
        Video(rr).h_avg.lin = [h_crest_lin_avg(xshift+1:end); NaN(xshift,1)];
        Video(rr).h_avg.nlin = [h_crest_nlin_avg(xshift+1:end); NaN(xshift,1); ];
        Video(rr).h_avg.bp = [h_crest_bp_avg(xshift+1:end); NaN(xshift,1); ];
    elseif isnan(xshift)
        Video(rr).h_avg.lin = h_crest_lin_avg;
        Video(rr).h_avg.nlin = h_crest_nlin_avg;
        Video(rr).h_avg.bp = h_crest_bp_avg;
    end

    %%% Mean gamma(x)
    gamma = zeros(length(Video(rr).x10),length(Video(rr).bp));
    for bb = 1:length(Video(rr).bp)
        if Video(rr).bp(bb) == 5001
            gamma(:,bb)=NaN;
        elseif isnan(Video(rr).bp(bb))
            gamma(:,bb)=NaN;
        else
            gamma(1:5001-Video(rr).bp(bb),bb) = 0.42;
        end
    end
    gamma(:,isnan(gamma(1,:)))=[];
    Video(rr).gamma = gamma;
    Video(rr).gamma_mean = mean(gamma,2);

    %%% Interpolate to 0.5m grid
    Video(rr).x5 = Video(rr).x10(1:5:end);
    Video(rr).h_avg.lin = interp1(Video(rr).x10, Video(rr).h_avg.lin, Video(rr).x5);
    Video(rr).h_avg.nlin = interp1(Video(rr).x10, Video(rr).h_avg.nlin, Video(rr).x5);
    Video(rr).h_avg.bp = interp1(Video(rr).x10, Video(rr).h_avg.bp, Video(rr).x5);
    
    %%% Cutoff offshore based on pixel resolution
    id = find([cutoff.date] == Video(rr).date & [cutoff.hover] == Video(rr).flight);
    id_mop = find(cutoff(id).mopgrid == round(Video(rr).mop));
    Video(rr).h_avg.lin(find(Video(rr).x5 == -cutoff(id).crest_track(id_mop)):end) = NaN;
    Video(rr).h_avg.nlin(find(Video(rr).x5 == -cutoff(id).crest_track(id_mop)):end) = NaN;
    Video(rr).h_avg.bp(find(Video(rr).x5 == -cutoff(id).crest_track(id_mop)):end) = NaN;

    %% Create composites
    
    [Video] = create_composite_bathys(Video, rr, gamma_H);
    
    %%% Variable surfzone RMSE
    lim = Video(rr).lims;
    
    % Inner Surfzone Stats
    [Video(rr).Error.RMSE_insz.cb, Video(rr).Error.Skill_insz.cb, Video(rr).Error.Bias_insz.cb] = calc_errors(Video(rr).survey.z, Video(rr).cbathy.z,[lim(1) lim(2)]);
    [Video(rr).Error.RMSE_insz.cb_hErr, Video(rr).Error.Skill_insz.cb_hErr, Video(rr).Error.Bias_insz.cb_hErr] = calc_errors(Video(rr).survey.z, Video(rr).cbathy.cbathy_hErr,[lim(1) lim(2)]);
    [Video(rr).Error.RMSE_insz.cb_gamma, Video(rr).Error.Skill_insz.cb_gamma, Video(rr).Error.Bias_insz.cb_gamma] = calc_errors(Video(rr).survey.z, Video(rr).cbathy.cbathy_gamma,[lim(1) lim(2)]);
    
    [Video(rr).Error.RMSE_insz.lin, Video(rr).Error.Skill_insz.lin, Video(rr).Error.Bias_insz.lin] = calc_errors(Video(rr).survey.z, interp1(Video(rr).x5, Video(rr).h_avg.lin, Video(rr).x10), [lim(1) lim(2)]);
    [Video(rr).Error.RMSE_insz.nlin, Video(rr).Error.Skill_insz.nlin, Video(rr).Error.Bias_insz.nlin] = calc_errors(Video(rr).survey.z, interp1(Video(rr).x5, Video(rr).h_avg.nlin, Video(rr).x10), [lim(1) lim(2)]);
    [Video(rr).Error.RMSE_insz.bp, Video(rr).Error.Skill_insz.bp, Video(rr).Error.Bias_insz.bp] = calc_errors(Video(rr).survey.z, interp1(Video(rr).x5, Video(rr).h_avg.bp, Video(rr).x10), [lim(1) lim(2)]);
    
    [Video(rr).Error.RMSE_insz.comp_hErr, Video(rr).Error.Skill_insz.comp_hErr, Video(rr).Error.Bias_insz.comp_hErr] = calc_errors(Video(rr).survey.z, Video(rr).composite.cbathy_hErr, [lim(1) lim(2)]);
    [Video(rr).Error.RMSE_insz.comp_gamma, Video(rr).Error.Skill_insz.comp_gamma, Video(rr).Error.Bias_insz.comp_gamma] = calc_errors(Video(rr).survey.z, Video(rr).composite.cbathy_gamma, [lim(1) lim(2)]);
    [Video(rr).Error.RMSE_insz.comp_nlin, Video(rr).Error.Skill_insz.comp_nlin, Video(rr).Error.Bias_insz.comp_nlin] = calc_errors(Video(rr).survey.z, Video(rr).composite.cbathy_nlin, [lim(1) lim(2)]);
    [Video(rr).Error.RMSE_insz.comp_CT, Video(rr).Error.Skill_insz.comp_bp, Video(rr).Error.Bias_insz.comp_bp] = calc_errors(Video(rr).survey.z, Video(rr).composite.cbathyCT, [lim(1) lim(2)]);
    
    % Breaking region Stats
    [Video(rr).Error.RMSE_break.cb, Video(rr).Error.Skill_break.cb, Video(rr).Error.Bias_break.cb] = calc_errors(Video(rr).survey.z, Video(rr).cbathy.z,[lim(2) lim(3)]);
    [Video(rr).Error.RMSE_break.cb_hErr, Video(rr).Error.Skill_break.cb_hErr, Video(rr).Error.Bias_break.cb_hErr] = calc_errors(Video(rr).survey.z, Video(rr).cbathy.cbathy_hErr,[lim(2) lim(3)]);
    [Video(rr).Error.RMSE_break.cb_gamma, Video(rr).Error.Skill_break.cb_gamma, Video(rr).Error.Bias_break.cb_gamma] = calc_errors(Video(rr).survey.z, Video(rr).cbathy.cbathy_gamma,[lim(2) lim(3)]);
    
    [Video(rr).Error.RMSE_break.lin, Video(rr).Error.Skill_break.lin, Video(rr).Error.Bias_break.lin] = calc_errors(Video(rr).survey.z, interp1(Video(rr).x5, Video(rr).h_avg.lin, Video(rr).x10), [lim(2) lim(3)]);
    [Video(rr).Error.RMSE_break.nlin, Video(rr).Error.Skill_break.nlin, Video(rr).Error.Bias_break.nlin] = calc_errors(Video(rr).survey.z, interp1(Video(rr).x5, Video(rr).h_avg.nlin, Video(rr).x10), [lim(2) lim(3)]);
    [Video(rr).Error.RMSE_break.bp, Video(rr).Error.Skill_break.bp, Video(rr).Error.Bias_break.bp] = calc_errors(Video(rr).survey.z, interp1(Video(rr).x5, Video(rr).h_avg.bp, Video(rr).x10), [lim(2) lim(3)]);
    
    [Video(rr).Error.RMSE_break.comp_hErr, Video(rr).Error.Skill_break.comp_hErr, Video(rr).Error.Bias_break.comp_hErr] = calc_errors(Video(rr).survey.z, Video(rr).composite.cbathy_hErr, [lim(2) lim(3)]);
    [Video(rr).Error.RMSE_break.comp_gamma, Video(rr).Error.Skill_break.comp_gamma, Video(rr).Error.Bias_break.comp_gamma] = calc_errors(Video(rr).survey.z, Video(rr).composite.cbathy_gamma, [lim(2) lim(3)]);
    [Video(rr).Error.RMSE_break.comp_nlin, Video(rr).Error.Skill_break.comp_nlin, Video(rr).Error.Bias_break.comp_nlin] = calc_errors(Video(rr).survey.z, Video(rr).composite.cbathy_nlin, [lim(2) lim(3)]);
    [Video(rr).Error.RMSE_break.comp_CT, Video(rr).Error.Skill_break.comp_bp, Video(rr).Error.Bias_break.comp_bp] = calc_errors(Video(rr).survey.z, Video(rr).composite.cbathyCT, [lim(2) lim(3)]);
    
    % Surfzone Stats
    [Video(rr).Error.RMSE_sz.cb, Video(rr).Error.Skill_sz.cb, Video(rr).Error.Bias_sz.cb] = calc_errors(Video(rr).survey.z, Video(rr).cbathy.z,[lim(1) lim(3)]);
    [Video(rr).Error.RMSE_sz.cb_hErr, Video(rr).Error.Skill_sz.cb_hErr, Video(rr).Error.Bias_sz.cb_hErr] = calc_errors(Video(rr).survey.z, Video(rr).cbathy.cbathy_hErr,[lim(1) lim(3)]);
    [Video(rr).Error.RMSE_sz.cb_gamma, Video(rr).Error.Skill_sz.cb_gamma, Video(rr).Error.Bias_sz.cb_gamma] = calc_errors(Video(rr).survey.z, Video(rr).cbathy.cbathy_gamma,[lim(1) lim(3)]);
    
    [Video(rr).Error.RMSE_sz.lin, Video(rr).Error.Skill_sz.lin, Video(rr).Error.Bias_sz.lin] = calc_errors(Video(rr).survey.z, interp1(Video(rr).x5, Video(rr).h_avg.lin, Video(rr).x10), [lim(1) lim(3)]);
    [Video(rr).Error.RMSE_sz.nlin, Video(rr).Error.Skill_sz.nlin, Video(rr).Error.Bias_sz.nlin] = calc_errors(Video(rr).survey.z, interp1(Video(rr).x5, Video(rr).h_avg.nlin, Video(rr).x10), [lim(1) lim(3)]);
    [Video(rr).Error.RMSE_sz.bp, Video(rr).Error.Skill_sz.bp, Video(rr).Error.Bias_sz.bp] = calc_errors(Video(rr).survey.z, interp1(Video(rr).x5, Video(rr).h_avg.bp, Video(rr).x10), [lim(1) lim(3)]);
    
    [Video(rr).Error.RMSE_sz.comp_hErr, Video(rr).Error.Skill_sz.comp_hErr, Video(rr).Error.Bias_sz.comp_hErr] = calc_errors(Video(rr).survey.z, Video(rr).composite.cbathy_hErr, [lim(1) lim(3)]);
    [Video(rr).Error.RMSE_sz.comp_gamma, Video(rr).Error.Skill_sz.comp_gamma, Video(rr).Error.Bias_sz.comp_gamma] = calc_errors(Video(rr).survey.z, Video(rr).composite.cbathy_gamma, [lim(1) lim(3)]);
    [Video(rr).Error.RMSE_sz.comp_nlin, Video(rr).Error.Skill_sz.comp_nlin, Video(rr).Error.Bias_sz.comp_nlin] = calc_errors(Video(rr).survey.z, Video(rr).composite.cbathy_nlin, [lim(1) lim(3)]);
    [Video(rr).Error.RMSE_sz.comp_CT, Video(rr).Error.Skill_sz.comp_bp, Video(rr).Error.Bias_sz.comp_bp] = calc_errors(Video(rr).survey.z, Video(rr).composite.cbathyCT, [lim(1) lim(3)]);
    
    % Full Profile Stats
    [Video(rr).Error.RMSE_full.comp_hErr, Video(rr).Error.Skill_full.comp_hErr, Video(rr).Error.Bias_full.comp_hErr] = calc_errors(Video(rr).survey.z, Video(rr).composite.cbathy_hErr, [1 5001]);
    [Video(rr).Error.RMSE_full.comp_gamma, Video(rr).Error.Skill_full.comp_gamma, Video(rr).Error.Bias_full.comp_gamma] = calc_errors(Video(rr).survey.z, Video(rr).composite.cbathy_gamma, [1 5001]);
    [Video(rr).Error.RMSE_full.comp_nlin, Video(rr).Error.Skill_full.comp_nlin, Video(rr).Error.Bias_full.comp_nlin] = calc_errors(Video(rr).survey.z, Video(rr).composite.cbathy_nlin, [1 5001]);
    [Video(rr).Error.RMSE_full.comp_CT, Video(rr).Error.Skill_full.comp_bp, Video(rr).Error.Bias_full.comp_bp] = calc_errors(Video(rr).survey.z, Video(rr).composite.cbathyCT, [1 5001]);
    
    % Offshore State
    [Video(rr).Error.RMSE_offshore.comp_hErr, Video(rr).Error.Skill_offshore.comp_hErr, Video(rr).Error.Bias_offshore.comp_hErr] = calc_errors(Video(rr).survey.z, Video(rr).composite.cbathy_hErr, [lim(3) 5001]);
    [Video(rr).Error.RMSE_offshore.comp_gamma, Video(rr).Error.Skill_offshore.comp_gamma, Video(rr).Error.Bias_offshore.comp_gamma] = calc_errors(Video(rr).survey.z, Video(rr).composite.cbathy_gamma, [lim(3) 5001]);
    [Video(rr).Error.RMSE_offshore.comp_nlin, Video(rr).Error.Skill_offshore.comp_nlin, Video(rr).Error.Bias_offshore.comp_nlin] = calc_errors(Video(rr).survey.z, Video(rr).composite.cbathy_nlin, [lim(3) 5001]);
    [Video(rr).Error.RMSE_offshore.comp_CT, Video(rr).Error.Skill_offshore.comp_bp, Video(rr).Error.Bias_offshore.comp_bp] = calc_errors(Video(rr).survey.z, Video(rr).composite.cbathyCT, [lim(3) 5001]);
    
end

% Save Data

Video_bathy = Video;
save(fullfile(data_dir, 'Video_bathy.mat'), 'Video_bathy', '-v7.3')
