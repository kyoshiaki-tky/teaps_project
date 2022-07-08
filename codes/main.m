function main
% Obtain Biologically stable and resilient parameter set for SBML model.
% REQUIREMENT: Parallel computing toolbox, Optimization toolbox,
%              SUNDIALS CVode solver toolbox(ver 2.6.0)
% In the process, a calculation module is added to path automatically,
% which is finally removed from path when successfully completed.

% RESULT: '../results/YEAR_MONTH/MODEL_NAME/RUN_date_time/Folder_com/'
%         "FinalX_fit.mat" includes finaly obtained parameter sets.
%         "modelFile.zip" is a zipped model file used in TEAPS calc. 
%         In 'RESULT_FOLDER/process/', files in calc process is recorded.

% Preparing a Model folder including Model XML file, Parameters.csv,
% Compartment.csv, odeSBML.m (ode description file), JacbSBML.m (jacobian descprition file)

% Throughout this script, state variables are written as "u"
% and kinetic parameters are written as "x"

currentDir = pwd;
Folder_com = ['new_code'];   % shortcomment, added on result folder name. 
model_name = ['NFKb5_2'];% specify model name in the "../Model/" folder.
                    % Cell Designer XML model should have
                    % same name in this folder.

%% Setting 
% initialization
Setting = struct;

%% Target condition setting

% describe a target fixed point in 'state_variables.csv' in the model folder.

% for calc O_basin,  target size of basin stability
Setting.target_basin_size = 0.1;  % considering the nature of the target system. our default is 0.1.

% for calc O_relax, max of non zero eigne valiue at the target fixed point 
Setting.target_max_eig = -0.3;  % aiming to observe convergence until time t=100

Setting.normEXP = true;     % using norm dU as exponential scale in objective function. default: true.
%Setting.normEXP = false;   % if using norm dU as normal scale, uncomment this line 

%% CNM Settings
Setting.Model = model_name;
Setting.nDivSamp = 2;  % "number of subcluster"; more than 2, default 10, 
Setting.samp= 400;     % recommendation: >=200, total sampling point number = nDivSamp x samp; 

%% Globalizaiton setting
% degree of tangential noise in globalized LBFGS
Setting.maxTan_shift = 10; %1~10 (arbitrary), default 10
% degree of noise in expanding of obtained sets
Setting.noise = 5;         %1~10 (arbitrary), default 5
% exponetial weight of O_fix in globalization process
Setting.weight_fix = 2;    % 2 (usually fix)
% exponetial weight of O_relax in globalization process
Setting.weight_relax = 1.2;  % 1~2 (starting with 2, if teaps returns no fit parameter sets, decreased this value.) 
% set exponetial weight of O_relax in globalization process
Setting.weight_basin = 1.2;  % usually set the same value with weight for O_relax

%% Serch space stringency
% 0: In the final optimization step, solution search starting from the
% search space, but allows solution outside serach space (recommended for large scale models)
% 1: try to find solution inside serach space in final optimization step
% anyway (using trust-regiondogleg method int the final solution search)
% Setting.SSS = 0;
Setting.SSS = 1;

%% Convergence criteria
% alpha: significant level in wilcoxon's rank sum test to test whether
% median values for all parameters are shifted by addition of parameter set
% obtained in latest iteration.
% beta: threshold of minimun inclusion rate of parameter values found in
% the previously found parameter distribution range for all parameters.
Setting.alpha = 0.1;
Setting.beta = 0.99;


%% other fixed setting
Setting.odeSolvers = {@CVode}   % ODE solver: fixed


%% Prepare folders
[Setting,calcFolder,resultFolder2,calcID] = prepareFolder(Setting);
disp(['calcualtion ID: ', calcID])
resultFolder = [resultFolder2, Folder_com,'/'];
mkdir(resultFolder);

zip('../modelFile.zip',[calcFolder '../Model/' Setting.Model]); 
movefile('../modelFile.zip',resultFolder);

%% Prepare parameters

cd('import');
parameters = ...
    str2double(importStr(['../../Model/',Setting.Model, '/','Parameters.csv'], ','))';  % importing Parameters.csv
cd('../');

Setting.resultFolder = [resultFolder, 'process/'];
mkdir(Setting.resultFolder)
csvwrite([Setting.resultFolder 'Parameters.csv'], parameters)
save([Setting.resultFolder 'Setting.csv'], 'Setting')
cd(currentDir)


%% Submit job
mainCal(Setting)

%% Complete
disp([char(10), 'Completed!',char(10)])
disp(['ResultFolder:', char(10), resultFolder])

cd(currentDir)

end

function [Setting,calcFolder,resultFolder,calcID] = prepareFolder(Setting)

calcFolder = [cd '/'];
cDate = datestr(now, 'yyyymmdd-HHMM');
cMonth = datestr(now, 'yyyymm');
calcID = cDate;
isFoldExist = 1;

while isFoldExist
    resultFolder = ...
    [pwd, '/../results/' cMonth '/' Setting.Model ...
    '/' calcID '/'];

    [~,~,messid] = mkdir(resultFolder);
    
    %Change folder name if the folder already exist
    if strcmp(messid,'MATLAB:MKDIR:DirectoryExists')==0
        isFoldExist = 0;
    else
        isFoldExist = isFoldExist + 1;
        calcID = [cDate '_' num2str(isFoldExist)];
    end
end

end
