%% Pixel resolution
% Using the CIRN Station Design toolbox
% 
%% Input 
%   file naming convention - YYYYMMDD_location_hover_MOP
%       if not, then need to change use of aa in 'Get stats from metadata'
%   Update camera specs based on given camera
%   Have metadata on Altitude, Yaw, Pitch and Roll of camera
%   Know grid you need to use for after (line 120)
%
%
%% Installation
%
%   requires toolbox from CIRN 
%       <https://github.com/Coastal-Imaging-Research-Network/Support-Routines>
%       <https://github.com/Coastal-Imaging-Research-Network/Station-Design-Toolbox>
%
%% Output 
%   cutoff (structure)
%       date
%       hover #
%       mopgrid - x-location
%       crest_track - for every MOP point, what is the offshore distance
%       cutoff based on pixel = 2*gridding for crest-tracking (0.2m)
%       cbathy - for every MOP point, what is the offshore distance
%       cutoff based on pixel = 2*gridding for cbathy (2m)
% 
% 
%%
function [cutoff] = pixel_res(files, data_dir, local_dir)
    %% Camera Specs:
    nameCam = 'THEIA';  % String that represents the name of the camera
                                 %   E.g. 'Flir Boson 640',...
    
    numCam = 1; % Number of cameras to be used.
    
    focalLength = 8.8;  % Focal length of the camera lens in mm.
    
    NU = 3840;   % Width of the chip (A.K.A. sensor) in pixels
    NV = 2160;   % Height of the chip in pixels
    
    hfov = 76.25;
    vfov = 47.64;
    
    roll=0;
    txyz=[0,0,0];
    %% Get stats from metadata
    clear stats aa
    for ff = length(files):-1:1
        aa(ff,:) = strsplit(files(ff).name, '_');
    end

    files = unique(append(aa(:,1), '_', aa(:,2), '_',aa(:,3)));



    for ff = 1:length(files)
        aa = strsplit(char(string(files(ff))), '_');
        cd([data_dir '/' char(string(aa(1))) '_' char(string(aa(2))) '/' char(string(aa(3)))])
        
        stats(ff).date = char(string(aa(1)));
        stats(ff).hover = str2double(char(string(aa(3))));
        C = readtable([char(string(aa(3))) '.csv']);
        % Elevation
        alt = C.GPSAltitude;
        if contains(aa{1,1}(3), '.')
            for ll = 1:size(alt,1)
                ab(ll)=str2double(alt{ll,1}(1:4));
            end
            stats(ff).elevation(nn) = max(ab);
        else
            for ll = 1:size(alt,1)
                ab(ll)=str2double(alt{ll,1}(1:2));
            end
            stats(ff).elevation = max(ab);
        end
        % Rotation
        stats(ff).heading = 360 + C.CameraYaw(end);
        stats(ff).roll = C.CameraRoll(end);
        stats(ff).tilt = 90 + C.CameraPitch(end);
    
    end
    cd(local_dir);
    %% Compute range
    for ff = 1:length(stats)
            [XYZVertices, res, dcRange, daRange, dcProj, daProj] = computeFootprint(hfov, vfov, stats(ff).heading, stats(ff).tilt, stats(ff).roll, NU, NV, [0,0, stats(ff).elevation], txyz);
            [x,y]=meshgrid(res{1,1}.x, res{1,1}.y);
            stats(ff).res.x=x+100;
            stats(ff).res.y=y;
            stats(ff).res.daRange = daRange;
            stats(ff).res.dcRange = dcRange;
    end
    
    % find range values for crest-tracking and cBathy
    for ff = 1:length(stats)
            stats(ff).res.range_ct = NaN(length(stats(ff).res.x), length(stats(ff).res.x));
            stats(ff).res.range_ct(find(stats(ff).res.daRange < 0.2 & stats(ff).res.dcRange < 0.2)) = 1;
    
            stats(ff).res.range_cb = NaN(length(stats(ff).res.x), length(stats(ff).res.x));
            stats(ff).res.range_cb(find(stats(ff).res.daRange < 2 & stats(ff).res.dcRange < 2)) = 1;
    end
    
    
    mopgrid = [-300:100:300]; 
    idy_ct = NaN(length(stats), length(mopgrid));
    idy_cb = NaN(length(stats), length(mopgrid));
    for ff = 1:length(stats)
            for mm = 1:length(mopgrid)
                mop=mopgrid(mm);
                idx = find(stats(ff).res.x(1,:)==mop);
                idys = stats(ff).res.y(find(stats(ff).res.range_ct(idx,:)==1),idx);
                if ~isempty(idys)
                    idy_ct(ff,mm)=min(idys);
                end
                idys = stats(ff).res.y(find(stats(ff).res.range_cb(idx,:)==1),idx);
                if ~isempty(idys)
                    idy_cb(ff,mm)=min(idys);
                end
            end
    end
    
    %% Put into cutoff variable
    clear cutoff
    jj=0;
    for ff = 1:length(stats)
            jj=jj+1;
            cutoff(jj).date = str2num(stats(ff).date(1:8));
            cutoff(jj).hover = stats(ff).hover;
            cutoff(jj).mopgrid = mopgrid/100+582; % grid around MOP 582 - Torrey Pines
            cutoff(jj).crest_track = squeeze(idy_ct(ff,:));
            cutoff(jj).cbathy = squeeze(idy_cb(ff,:));
    end
    save([data_dir 'cutoff_pixres.mat'], 'cutoff')
end