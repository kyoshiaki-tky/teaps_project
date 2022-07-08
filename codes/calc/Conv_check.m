function Teaps_output = Conv_check( Teaps_output, Model, u_info,resultFolder)
%   checking convergence 

homeDir = pwd;

%% initialization
uconst = u_info.uconst;
uinitial=u_info.uinitial;
% uinitial=cellfun(@str2num,u_info.uinitial{1});
uConv_idx = find(u_info.uConv);

%% Calculation
cd('../calc');

Xvalue = Teaps_output.finalX;
nPS = size(Xvalue,2);

ulength = length(uConv_idx);

% if strcmp(Model, 'T8')
%     ulength = ulength-2;
% elseif (strcmp(Model, 'SS5')||strcmp(Model, 'SS6'))
% else
%     ulength = ulength-1;
% end 
    
dU_norm = zeros(nPS,1);
NZLM = zeros(nPS,1);

for kc2=1:nPS
    try
        x_1set = Xvalue(:,kc2);
        [dU_norm(kc2) NZLM(kc2)] = Conv_calc( uinitial, x_1set, Model, uconst); %#ok<NCOMMA>
    catch
        dU_norm(kc2) = nan;
        NZLM(kc2) = nan;
    end
end


% check fitness: set criteria
fit_norm = (dU_norm < 1e-7);
fit_NZLM = (NZLM < 0);

Teaps_output.fit_idx = (fit_norm & fit_NZLM);

save([resultFolder,'Conv.mat'])
cd(homeDir)
end
