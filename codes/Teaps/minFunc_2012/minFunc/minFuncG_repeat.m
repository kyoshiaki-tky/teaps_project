function [x,f,exitflag,output] = minFuncG_repeat(funObj,funObj2, x0,options,tanMaxLength,repN,ext,x_lowlim,x_uplim,x_ori2,varargin)
% repeating minFunG_0 for repN times.

x0_2 = x_ori2(:,(sum(isinf(x_ori2))==0));
X_mean = (mean(x0_2,2));
% X_median = (median(x0_2,2));

for k= 1:repN
    if k~=1
        x0_ori = x0;
%         disp([num2str(k), 'before  ext vector formation']);
        extVector = (x-x0_ori)/ norm(x-x0_ori);
        if sum(isnan(extVector)) || sum(isinf(extVector))
            extVector = rand(size(extVector))-0.5;
            extVector = extVector/norm(extVector);
        end
%         disp([num2str(k), 'after  ext vector formation']);
        
        axis_idx = 1;
        while abs(extVector(axis_idx)) < 1e-6
            axis_idx = axis_idx+1;
        end
%         disp([num2str(k), 'after  idx selection']);        
        idx_logical = zeros(size(extVector));
        idx_logical(axis_idx) = true; 
        extTanVector = rand(size(extVector));
        extTanVector(axis_idx) = -dot(extTanVector(~idx_logical),extVector(~idx_logical))/extVector(axis_idx);
        extTanVector = extTanVector./norm(extTanVector);   % normalized tangential vector

        width = log(x_uplim)-log(x_lowlim);
        width(width>20) = 20; 

        x0 = x + ext.*rand(1).*(rand(1).*extVector + extTanVector).*exp(-((x-X_mean)./(width/2)).^2);


    end

    % limitation of serach space
    % for nan
    x0(isnan(x0)) = rand(size(x0(isnan(x0)))).*(log(x_uplim(isnan(x0)))-log(x_lowlim(isnan(x0)))) + log(x_lowlim(isnan(x0)));
    % for x< low_lim
    x0(x0<log(x_lowlim)) = rand(size( x0(x0<log(x_lowlim)) )).*(log(x_uplim(x0<log(x_lowlim)))-log(x_lowlim(x0<log(x_lowlim))))/2 + (log(x_uplim(x0<log(x_lowlim)))+log(x_lowlim(x0<log(x_lowlim))))/2;
    % for x< up_lim
    x0(x0>log(x_uplim)) = -rand(size( x0(x0>log(x_uplim)) )).*(log(x_uplim(x0>log(x_uplim)))-log(x_lowlim(x0>log(x_uplim))))/2 + (log(x_uplim(x0>log(x_uplim)))+log(x_lowlim(x0>log(x_uplim))))/2;    

    % possible to implement two onjective dunction for automatic error
    % recovery. if an error occurs, returns original value. 
    try
%        'obj1'
        if isempty(varargin)
            [x,f,exitflag,output] = minFunc_global(funObj,x0,options,tanMaxLength,x_lowlim,x_uplim);
        else
            [x,f,exitflag,output] = minFunc_global(funObj,x0,options,tanMaxLength,x_lowlim,x_uplim,varargin{:});
        end
    catch
        if ~(isempty(funObj2))
            try
    %            'obj2'
                if isempty(varargin)
                    [x,f,exitflag,output] = minFunc_global(funObj2,x0,options,tanMaxLength,x_lowlim,x_uplim);
                else
                    [x,f,exitflag,output] = minFunc_global(funObj2,x0,options,tanMaxLength,x_lowlim,x_uplim,varargin{:});
                end
            catch
    %            'noshift'
                x = x0_ori;
            end
        else
            x = x0_ori;            
        end
    end
end

end