%% Create Composite Profiles  
% from cBathy and crest-tracking surfzone bathymetry
% 
%% Input 
%   Video (structure) - from bathy_inversion.m
%   rr - index of which value in Video structure to use
% 
%% Output 
%   Video with composite bathymetries 
%       composite.cbathy_hErr - cBathy with hErr > 0.5m removed and interpolated over
%       composite.cbathy_gamma - cBathy with gamma-derived breaking region removed and interpolated over
%       composite.cbathy_nlin - cBathy_gamma with nonlinear correction gamma(x) applied
%       composite.cbathyCT - BP crest-tracking in the surfzone and cBathy offshore
% 
%% Copyright 
% Athina Lange 2022
%
%%
function [Video] = create_composite_bathys(Video, rr, gamma)
    window = (50/2)/0.1;
%% Remove any cBathy hErr > 0.5 and interp
    Video(rr).cbathy.cbathy_hErr = Video(rr).cbathy.z;
    Video(rr).cbathy.cbathy_hErr(Video(rr).cbathy.zerr > 0.5) = NaN;
    x10_temp=Video(rr).x10;
    x10_temp(isnan(Video(rr).cbathy.cbathy_hErr))=[];Video(rr).cbathy.cbathy_hErr(isnan(Video(rr).cbathy.cbathy_hErr))=[];
    if ~isempty(x10_temp)
        Video(rr).cbathy.cbathy_hErr = interp1(x10_temp, Video(rr).cbathy.cbathy_hErr, Video(rr).x10);
    else
        Video(rr).cbathy.cbathy_hErr = NaN(size(Video(rr).x10,1), size(Video(rr).x10,2));
    end
%% Get cBathy with gamma criteron - including offshore and onshore breaking cutoff
    gamma_profile = Video(rr).gamma_mean; % Not applying xshift because was computed from unshift data - cBathy also unshifted
    aa = find(min(abs(gamma_profile-(gamma-0.01)))==abs(gamma_profile-(gamma-0.01))); id_shore = aa(end); % onshore
    aa = find(round(gamma_profile,2)==0); 
    if isempty(aa)
        aa = find(round(gamma_profile,2)==min(round(gamma_profile,2)));
    end
    
    id_ocean = aa(1)+window; % offshore + window/2
     
    cbathy_clean = Video(rr).cbathy.cbathy_hErr; xcb_temp = Video(rr).x10'; cbathy_clean(id_shore+1:id_ocean)=NaN; 
    xcb_temp(isnan(cbathy_clean))=[]; cbathy_clean(isnan(cbathy_clean))=[];
    if ~isempty(xcb_temp)
        Video(rr).cbathy.cbathy_gamma = interp1(xcb_temp, cbathy_clean, Video(rr).x10);
    else
        Video(rr).cbathy.cbathy_gamma = NaN(size(Video(rr).x10,1), size(Video(rr).x10,2));
    end
%% Foreshore Beach
    date = Video(rr).date;
    tide = Video(rr).min_tide;
    x = Video(rr).x10';
    z = Video(rr).survey.z';
    [~,id]=min(abs(z-tide)); % find when bathy reaches low tide line
    x(id+1:end)=[]; % only foreshore
    z(id+1:end)=[];
    x_upper = x'; z_upper = z';
    x_upper(isnan(z_upper))=[];z_upper(isnan(z_upper))=[];

%% cBathy Composite - only hErr < 0.5m
    x10_temp = Video(rr).x10;
    cb = Video(rr).cbathy.cbathy_hErr;
    x10_temp(isnan(cb))=[]; cb(isnan(cb))=[];

    if ~isempty(x10_temp)
        % if overlap between subaerial and cBathy -> prioritize subaerial survey
        if x10_temp(1) < x_upper(end)
            ii=find(min(abs(x10_temp - x_upper(end))) == abs(x10_temp - x_upper(end)));
            x10_temp(1:ii)=[]; cb(1:ii)=[];
        elseif x10_temp(1) == x_upper(end)
            x10_temp(1)=[]; cb(1)=[];
        end
        % combine subaerial survey and cBathy
        x = [x_upper; x10_temp];
        z = [z_upper; cb];      
        Video(rr).composite.cbathy_hErr = smoothdata(interp1(x, z, Video(rr).x10),'Gaussian',window);
        runupline = InterX([Video(rr).x10';Video(rr).composite.cbathy_hErr'],[x_upper';z_upper']); % correct for smoothing issue; find intersection point between upper beach and smoothed curve
        [~,i1]=min(abs(x_upper-runupline(1,1)));
        [~,i2]=min(abs(Video(rr).x10-runupline(1,1)));
        x = [x_upper(1:i1); Video(rr).x10(i2+1:end)];
        z = [z_upper(1:i1); Video(rr).composite.cbathy_hErr(i2+1:end)];
        Video(rr).composite.cbathy_hErr = interp1(x, z, Video(rr).x10);
        Video(rr).stitch.cbathy_hErr = x_upper(i1);
    else
        Video(rr).composite.cbathy_hErr = NaN(size(Video(rr).x10,1), size(Video(rr).x10,2));
    end
%% cBathy Composite - gamma
    x10_temp = Video(rr).x10;
    cb = Video(rr).cbathy.cbathy_gamma;
    x10_temp(isnan(cb))=[]; cb(isnan(cb))=[];

    if ~isempty(x10_temp)
        % if overlap between subaerial and cBathy -> prioritize subaerial survey
        if x10_temp(1) < x_upper(end)
            [~,ii]=min(abs(x10_temp - x_upper(end)));
            x10_temp(1:ii)=[]; cb(1:ii)=[];
        end
        % combine subaerial survey and cBathy
        x = [x_upper; x10_temp];
        z = [z_upper; cb];
                
        Video(rr).composite.cbathy_gamma = smoothdata(interp1(x, z, Video(rr).x10),'Gaussian',window);
        runupline = InterX([Video(rr).x10';Video(rr).composite.cbathy_gamma'],[x_upper';z_upper']); % correct for smoothing issue; find intersection point between upper beach and smoothed curve
        [~,i1]=min(abs(x_upper-runupline(1,1)));
        [~,i2]=min(abs(Video(rr).x10-runupline(1,1)));
        x = [x_upper(1:i1); Video(rr).x10(i2+1:end)];
        z = [z_upper(1:i1); Video(rr).composite.cbathy_gamma(i2+1:end)];
        Video(rr).composite.cbathy_gamma = interp1(x, z, Video(rr).x10);
        Video(rr).stitch.cbathy_gamma = x_upper(i1);
    else
        Video(rr).composite.cbathy_gamma = NaN(size(Video(rr).x10,1), size(Video(rr).x10,2));
    end
%% cBathy Composite - nonlinear
    x10_temp = Video(rr).x10;
    cb = Video(rr).cbathy.cbathy_gamma./ (1+Video(rr).gamma_mean);
    x10_temp(isnan(cb))=[]; cb(isnan(cb))=[];

     if ~isempty(x10_temp)
        % if overlap between subaerial and cBathy -> prioritize subaerial survey
        if x10_temp(1) < x_upper(end)
            [~,ii]=min(abs(x10_temp - x_upper(end)));
            x10_temp(1:ii)=[]; cb(1:ii)=[];
        end
        % combine subaerial survey and cBathy
        x = [x_upper; x10_temp];
        z = [z_upper; cb];
                
        Video(rr).composite.cbathy_nlin = smoothdata(interp1(x, z, Video(rr).x10),'Gaussian',window);
        runupline = InterX([Video(rr).x10';Video(rr).composite.cbathy_nlin'],[x_upper';z_upper']); % correct for smoothing issue; find intersection point between upper beach and smoothed curve
        [~,i1]=min(abs(x_upper-runupline(1,1)));
        [~,i2]=min(abs(Video(rr).x10-runupline(1,1)));
        x = [x_upper(1:i1); Video(rr).x10(i2+1:end)];
        z = [z_upper(1:i1); Video(rr).composite.cbathy_nlin(i2+1:end)];
        Video(rr).composite.cbathy_nlin = interp1(x, z, Video(rr).x10);
        Video(rr).stitch.cbathy_nlin = x_upper(i1);
    else
        Video(rr).composite.cbathy_nlin = NaN(size(Video(rr).x10,1), size(Video(rr).x10,2));
    end
%% BP Composite
    id_sz = 0; id_cb = 0;
    % Surfzone bathy
    x_sz = Video(rr).x10; z_sz = smoothdata(interp1(Video(rr).x5, Video(rr).h_avg.bp, Video(rr).x10), 'Gaussian', window);
    z_sz(1:find(max(x_upper)==x_sz)) = NaN; % find how far offshore subaerial survey region is valid
    z_sz(id_ocean:end) = NaN; % remove anything further offshore then beginning of breaking + cBathy buffer
    x_sz(isnan(z_sz)) = []; z_sz(isnan(z_sz)) = [];
    if ~isempty(x_sz)
        shore_id = x_sz(1);
        % if overlap between subaerial and crest-tracking -> prioritize subaerial survey
        if x_sz(1) < x_upper(end)
            [~,ii]=min(abs(x_sz - x_upper(end)));
            x_sz(1:ii)=[]; z_sz(1:ii)=[];
        end
        id_sz = 1;
    else
        shore_id = id_shore/10;
    end

    % Offshore bathy
    x_cb = Video(rr).x10;
    z_cb = Video(rr).cbathy.cbathy_hErr; z_cb(1:id_ocean) = NaN; 
    x_cb(isnan(z_cb)) = []; z_cb(isnan(z_cb)) = [];
    if ~isempty(x_sz); x_cut = x_sz(end); else; x_cut = x_upper(end);end
    if ~isempty(x_cb)
        if x_cb(1) < x_cut
            [~,ii]=min(abs(x_cb - x_cut));
            x_cb(1:ii)=[]; z_cb(1:ii)=[];
        end
        id_cb = 1;
    end
  
    % combine subaerial survey, BP crest-tacking, and cBathy
    x = [x_upper; x_sz; x_cb]; 
    z = [z_upper; z_sz; z_cb];

    Video(rr).composite.cbathyCT = smoothdata(interp1(x, z, Video(rr).x10),'Gaussian', window);
    runupline = InterX([Video(rr).x10';Video(rr).composite.cbathyCT'],[x_upper';z_upper']); % correct for smoothing issue; find intersection point between upper beach and smoothed curve
    [~,i1]=min(abs(x_upper-runupline(1,1)));
    [~,i2]=min(abs(Video(rr).x10-runupline(1,1)));
    x = [x_upper(1:i1); Video(rr).x10(i2+1:end)];
    z = [z_upper(1:i1); Video(rr).composite.cbathyCT(i2+1:end)];
    Video(rr).composite.cbathyCT = interp1(x, z, Video(rr).x10);
    Video(rr).stitch.cbathyCT_beach = x_upper(i1);
   
    if ~isempty(x_cb) &&  ~isnan(x_cb(1))
        Video(rr).stitch.cbathyCT_offshore = x_cb(1);
    elseif ~isempty(x_sz) &&  ~isnan(x_sz(end))
        Video(rr).stitch.cbathyCT_offshore = x_sz(end);
    end

    % index of [1st BP valid onshore point, onshore cutoff of breaking, offshore cutoff of breaking]
    Video(rr).lims = (round([shore_id*2; id_shore; id_ocean])); 
    if id_sz == 1; Video(rr).check.CT = 'Yes'; else Video(rr).check.CT = 'No';end
    if id_cb == 1; Video(rr).check.cBathy = 'Yes'; else Video(rr).check.cBathy = 'No';end

end
