function mainCal(Setting) %#ok<INUSD>

homeDir = pwd;

%% Setting

Model = Setting.Model;
resultFolder = Setting.resultFolder;
tic;

%% Prepare parameters
Numbers = struct;
cd([homeDir '/import']);
[Numbers InitialXY ParameterIndex u_info] = importParameters(Model, Numbers); %#ok<NCOMMA>
cd(homeDir);

Numbers.samp = Setting.samp;   
Numbers.samp = Numbers.samp*(Setting.nDivSamp);     % total sample number
Numbers.nDivSamp = Setting.nDivSamp;     % total sub-cluster number
Numbers.target_basin_size = Setting.target_basin_size;
Numbers.target_max_eig = Setting.target_max_eig;
Numbers.normEXP = Setting.normEXP;
Numbers.maxTan_shift = Setting.maxTan_shift;
Numbers.noise = Setting.noise;
Numbers.SSS = Setting.SSS;   % serach space stringecy
Numbers.weight_fix = Setting.weight_fix;
Numbers.weight_relax = Setting.weight_relax;
Numbers.weight_basin = Setting.weight_basin;
Numbers.alpha = Setting.alpha;
Numbers.beta = Setting.beta;

mkdir([resultFolder 'Numbers']);
save([resultFolder 'Numbers/Numbers.mat'], 'Numbers')


%% TEAPS Calculation

cd([homeDir '/Teaps']);
k=1;
Teaps_output(k) = {Teaps_run(InitialXY, Numbers, ParameterIndex, Model, u_info, resultFolder)};

% checking convergence of distribution
dist_conv = false;
alpha = Numbers.alpha;
beta = Numbers.beta;
while ~dist_conv
    k = k+1;
    Teaps_output(k) = {Teaps_run(InitialXY, Numbers, ParameterIndex, Model, u_info, resultFolder)};
    % checking median shift
    [dist_conv,Teaps_output_all] = dist_conv_judge(Teaps_output,alpha);
    % checking inclusion rate
    if dist_conv
        dist_conv = dist_conv_judge2(Teaps_output,beta);
    end
    if ~dist_conv
        disp('Distribution shifted by last itelation of TEAPS.') 
        disp(['next iteration number is #', num2str(k+1), '.'])
    else
        disp(['No significant distribution shift by iteration #', num2str(k), '.'])
        disp('Finished TEAPS iterations.')
    end
end



%% testing du_norm and NZLM
cd([homeDir '/calc']);
Teaps_output_all= Conv_check( Teaps_output_all, Model, u_info, resultFolder);
disp(['fitted parameter set: ', num2str(sum(Teaps_output_all.fit_idx)), ' sets (out of ', num2str(Numbers.samp*k), ' sets)'])

FinalX_fit = Teaps_output_all.finalX(:,Teaps_output_all.fit_idx);
save([resultFolder, '../FinalX_fit.mat'],'FinalX_fit')
save([resultFolder, 'Teaps_output.mat'],'Teaps_output')

%% simulation of fitted parameter set
% simulate with n_Sim parameter sets at maximum with TEAPS parameter sets randomly selected from
Tmax = 100;   % Simulate until Tmax
n_Tdiv = 200; % number of time points in figure plotting
n_Sim = 100;  % maximum number of parameter sets to be simulated
Rth = 0.1;    % max perturbation noise of initial state

sim_fitTC( Teaps_output_all, Model,u_info, Tmax, n_Tdiv , Rth, n_Sim, resultFolder);

%% Close
csvwrite([resultFolder 'tocf.csv'], tocf(6))
tocf(6)
cd(homeDir);
end
