%% Error calculation
%
%% Input
%   obs - ground-truth data that prediction is being compared to ([1 x n])
%   pred - prediction ([1 x n])
%   lim - index of limits of area that you want to compute statistics for ([id1 id2])
%
%% Output
%   r - RMSE
%   s - Skill
%   b - Bias
%
%% Copyright
% Athina Lange 2022
%
%%
function [r,s,b] = calc_errors(obs, pred, lim)
    tempx = pred; tempx = tempx(lim(1):lim(2));
    ax = obs; ax = ax(lim(1):lim(2));
    
    ax(find(isnan(tempx)))=[];
    tempx(find(isnan(tempx)))=[];
    tempx(find(isnan(ax)))=[];
    ax(find(isnan(ax)))=[];
    [r,s,b]=error(ax, tempx);

    function [rmse, skill, bias] = error(obs, pred)
        rmse = sqrt(sum((pred - obs).^2)/length(obs));
        skill = 1 - sum((pred - obs).^2)./sum((obs - mean(obs)).^2);
        bias = sum(pred - obs) / length(obs);
    end

end