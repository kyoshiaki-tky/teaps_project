function Teaps_output = Teaps_run(InitialXY, Numbers, ParameterIndex, Model, u_info, resultFolder)
% parameter distribution determined by CNM. Optimized by fixed point, jacobian
% also restriction by lamda value (in order to guarantee the speed of
% relaxation)

format shorte

curDir = pwd;
cd(['../../Model/', Model])
odeF2= @odeSBML;
JacbF2 = @JacbSBML;
MB_constraint2= @MB_constraint;
cd(curDir)

% define calculation tolerance 
cal_tol = 1e-8;    

nDivSamp = Numbers.nDivSamp;
minimum_value(:,1:Numbers.xfree) = InitialXY.xfree_import(1,:);
maximum_value(:,1:Numbers.xfree) = InitialXY.xfree_import(2,:);
interval = maximum_value-minimum_value;   % calc width of search space for each parameter (x)

% random generator initiation 
rng(now*1E3,'twister'); 
rand_val = rand(Numbers.xfree,Numbers.samp);
    

%% norm exponential    % chose objective function shape for norm of dU
normEXP = Numbers.normEXP;     % using norm dU as exponential scale in objective function


%% Target setting

% target fixed point
uconst = u_info.uconst;
uinitial=u_info.uinitial;

y_star_set(1,1) = 1;  % fixed point : answer for O_Fix (will be defined later)
y_star_set(2,1) = 1;  % basin       : answer for O_basin (will be defined later)
y_star_set(3,1) = 1;  % relax       : answer for O_relax (will be defined later)

NZLM_target = Numbers.target_max_eig;
Numbers_goal= 3;

%% Prepare parameters

x_lowlim = InitialXY.xfree_ss_import(1,:)';
x_uplim = InitialXY.xfree_ss_import(2,:)';

InitialXY.xfree_distribution = zeros(Numbers.xfree,Numbers.samp);
xmat = zeros(Numbers.xfree,Numbers.samp);
xmat(1:Numbers.xfree,:) = repmat(minimum_value',1,Numbers.samp)+(repmat(interval',1,Numbers.samp).*rand_val); 
InitialXY.X0 = xmat; 

InitialXY.xfree_distribution(1:Numbers.xfree,:) = xmat; 
X = zeros(Numbers.xfree,Numbers.samp,1);
X(:,:,1) = xmat;


Y = zeros(Numbers_goal,Numbers.samp,1);    % check later
A = zeros(Numbers_goal, Numbers.xfree, 1,nDivSamp);  % check later
l_not_valid = zeros(3,1); l_valid = zeros(3,1);


%% Prepare X and Y for FIRST calculation

% convert to normal scale value when exponential search
XforCalculation = prepareXforCalculation(Numbers,InitialXY,ParameterIndex,X(:,:,1));

% just display XforCalc top 5 rows.
XforCalculation(1:min(Numbers.samp,5),ParameterIndex.free)' %#ok<NOPRT>


%% Finding out the stable initial x sets
uNOTconst = ~uconst;
uNOTconstNum = sum(uNOTconst);

%% Addition of path to optimization module
addpath(genpath('minFunc_2012'));

%% MAIN ROUTINE
parpoolNum = maxNumCompThreads*2-1;
if parpoolNum > 20
    parpoolNum = parpoolNum-1;
end
try
    parpool(parpoolNum)
catch
    delete(gcp)
    parpool(parpoolNum)
end
cd(['../../Model/', Model])
addAttachedFiles(gcp, {'odeSBML.m','JacbSBML.m'})
cd(curDir)

N_obsP = 20;       % number of points of the obseravation point % for samall model 500;
x_store=[];           % initialization of parameer sets fitted to the collection criteria

%% stable region expansion setting
maxBasinSize = Numbers.target_basin_size;   % maximum perturbation
gradN=3;                                    % number of steps in expanding basin
fitsum0_count = 0;  % the count of sequential subroutines in which the target parameter sets have not been collected

for k_SS = 1:2:(2*(gradN-1)+1) % index for stable region size
    Y = zeros(Numbers_goal,Numbers.samp,1);   %initialization of Y
    A = zeros(Numbers_goal, Numbers.xfree, 1,nDivSamp);  %initialization of A
    x_store  = [];        % initialization of xfit in a routine with a certain size of observation are (initalized when ii is increased)
    u_fix = uinitial;  % importing target fixed point state, vector: parameter number x 1.
    u2    = repmat(u_fix,1,N_obsP);    % prep for making obsP as matrix
    u2ori = u2;                 % storing initial value matrix (original) 
       
   % allocating the memory space for follows (max 15000 iterations in CNM)
    duNorm = zeros(15000,4);
    MRJ    = zeros(15000,4);
    FPJ    = zeros(15000,4);

    k_CNM = 1;     % initialization idex for CNM subroutine
    if k_SS == 1
        disp('CNM process')
    end

    u2(uNOTconst,1:N_obsP) = repmat(u_fix(uNOTconst),1,N_obsP).* ...
        (1+(rand(length(u_fix(uNOTconst)),N_obsP)-0.5)/0.5*maxBasinSize*((k_SS)/(2*(gradN-1)+1))^2);
    u2 = u_constraint(u2, MB_constraint2);

    while (k_CNM<15001) && k_SS==1
        % observation point
        u2(uNOTconst,1:N_obsP) = repmat(u_fix(uNOTconst),1,N_obsP).* ...
            (1+(rand(length(u_fix(uNOTconst)),N_obsP)-0.5)/0.5*maxBasinSize*((k_SS)/(2*(gradN-1)+1))^2);
        u2 = u_constraint(u2, MB_constraint2);

        oldID = zeros(nDivSamp,Numbers.samp/(nDivSamp));
        X_k_subCL = zeros(Numbers.xfree, Numbers.samp/nDivSamp ,nDivSamp);
        Y_k_subCL = zeros(Numbers_goal, Numbers.samp/nDivSamp ,nDivSamp);
        fitIndex_subCL = zeros(nDivSamp, Numbers.samp/nDivSamp);

        parfor k1=1:Numbers.samp %#ok<*PFUNK>  %parfor
            x = XforCalculation(k1,:);
            
            % 1st indicator: Calculation of norm of time derivertive of state u (dU)  
            % norm of time derivertive of state u_fix 
            dU = odeF2(x,u_fix);            
            if normEXP   % exponential cal on norm od dU in obj function
                O_fix(1,k1) = exp(norm(dU(uNOTconst))/sqrt(uNOTconstNum));
                  % by dividing sqrt(number of U), the size of dU norm is 
                  % normalized to calcel out the effect of varibale number 
                  % uNOTconst: the index vector filtering not constant
                  % vatibale
            else         % normal scale calc on norm dU in obj fun 
                O_fix(1,k1) = norm(dU(uNOTconst))/sqrt(uNOTconstNum)+1; % for preventing division by 0, +1 addition
            end
            
            % 2nd indicator: NON-ZERO Eigenvalue of target steady state 
            % calculate  Jacobian matrix at state U using parameter set x.
            Jori = JacbF2(x, u_fix);
            try
                % calculation of Real eigenvalue of Jacobian matrix at
                % target steady state
                rj_ss = real(eig(Jori(uNOTconst, uNOTconst)));
                % ReLU(Real eig of Jacb matrix)+1
                O_relax(1,k1) = max([(max(rj_ss(abs(rj_ss)> cal_tol))-(NZLM_target)),0])+1;% for preventing division by 0, +1 addition
            catch
                % if jacb calc fails: return inf to exclude current 
                % parameter set in the further calculations 
                O_relax(1,k1) = Inf;
            end
            
            % 3rd indicator: Max Eigenvalue of focused stabale reagion
            % initialize max eigen value as -inf:
            O_basin(1,k1)=-Inf;
            % then, update maximum eigenvalue by iteratively calculating
            % maximum eigenvalue for each observation point in U2
            for ks=1:N_obsP
                u3 = u2(:,ks); %#ok<PFBNS>   % select one obsertvation point as u3
                % calculation of Jacobian matrix at state u3 using parameter set x 
                J = JacbF2( x, u3);
                try
                    % excluding rows and columns for constant variable
                    J_notConst = J(uNOTconst, uNOTconst);
                    % Calc real(eig(J+J_inv))
                    rj_u3_ks = real(eig( J_notConst + J_notConst' ));
                    % Calc current max of real eigenvalue of (J+J_inv)
                    mrj_current = max(rj_u3_ks);
                catch
                    % if jacb calc fails: return inf to exclude current 
                    % parameter set in the further calculations 
                    mrj_current = Inf;
                end
                % Comparison with recorded maxium vs current value
                % If current value is greater, then renew.
                O_basin(1,k1) = max(mrj_current,O_basin(1,k1));
            end
            
            O_basin(1,k1) = (O_basin(1,k1)>cal_tol).*O_basin(1,k1) +1;  % for preventing division by 0, +1 addition
        end
        %        Y_k_half =[fixed point; stability ; relaxation speed];        
        Y_k_CL =[O_fix; O_basin; O_relax];

        duNorm(k_CNM,:)=[max(log(O_fix)),median(log(O_fix)),mean(log(O_fix)),min(log(O_fix))];
        MRJ(k_CNM,:)=[log(abs(max(O_basin))),log(abs(median(O_basin))),log(abs(mean(O_basin))),log(abs(min(O_basin)))];
        FPJ(k_CNM,:)=[log(abs(max(O_relax))),log(abs(median(O_relax))),log(abs(mean(O_relax))),log(abs(min(O_relax)))];

        valid_idx2=find(sum(Y_k_CL~=Inf)); not_valid_idx2=find(sum(Y_k_CL==Inf));
        l_not_valid(1)=size(not_valid_idx2,2); l_valid(1)=size(valid_idx2,2);
        if l_valid(1)==0, 'error: max eigen value of all samples reached Inf', break,end %#ok<NOPRT>
        
        % Find NaN samples
        nan_index = find(isnan(sum(Y_k_CL,1)));
        l_nan=size(nan_index,2);
        if l_nan~=0, l_nan,error('y_out reached NaN'), end %#ok<NOPRT>
        Y(:,valid_idx2,1)=Y_k_CL(1:Numbers_goal,valid_idx2);

        if k_CNM ==1
            target_fix   = 1 + 1e-5;  % arbitrary set near the goal
            target_basin = 1 + 1e-5;  % arbitrary set near the goal
            target_relax = 1 + 1e-5;  % arbitrary set near the goal
            targetVal =[target_fix, target_basin, target_relax];
        end
        
        % consider renaming
        fitIndex1 = Y_k_CL(1,:) < targetVal(1);
        fitIndex2 = Y_k_CL(2,:) < targetVal(2);
        fitIndex3 = Y_k_CL(3,:) < targetVal(3);
        
        fitIndex = fitIndex1 & fitIndex2 & fitIndex3;
        xfit_current_iter = X(:,fitIndex,1);
        x_store = [x_store xfit_current_iter]; % storing parameter sets near target

        progress_bar_w = 50; % definition of progress bar width
        iterN = sprintf('iterN:%5.0f', k_CNM);
        if size(x_store,2)>=Numbers.samp
            if k_SS==1 && k_CNM< 20  % to find wider parameter set, CNM loop continues at least 20 iterations.
                percent = 100;
                perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
                message = [ iterN, char(10), perc,'[', repmat('=', 1, round(percent*progress_bar_w/100)), '=', repmat(' ', 1, progress_bar_w - round(percent*progress_bar_w/100)), ']'];
                rm_message = repmat(char(8),1,length(message)+1);
                disp([rm_message, message]);            
            else
                percent = 100;
                perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
                message = [ iterN,char(10), perc,'[', repmat('=', 1, round(percent*progress_bar_w/100)), '=', repmat(' ', 1, progress_bar_w - round(percent*progress_bar_w/100)), ']'];
                rm_message = repmat(char(8),1,length(message)+1);
                disp([rm_message, message]);

                [randNum randIndex] =sort(rand(1,size(x_store,2))); %#ok<NCOMMA,ASGLU>
                x_store_final = x_store(:,randIndex(1:Numbers.samp));
                X(:,:,1) = x_store_final;
                Teaps_output_raw(:,:,2*k_SS-1) = X(:,:,1); %#ok<AGROW>
                break
            end
        else
            percent = size(x_store,2)/Numbers.samp;
            perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
            message = [iterN, char(10), perc,'[', repmat('=', 1, round(percent*progress_bar_w/100)), '>', repmat(' ', 1, progress_bar_w - round(percent*progress_bar_w/100)), ']'];
            rm_message = repmat(char(8),1,length(message)+1);
            if k_CNM>1
                disp([rm_message, message]);
            else
                disp([char(10),message]);
            end            
        end

        %% Inverse transformation        
        % dividing to subclusters
        if nDivSamp>1
            inf_indexC = cell(1,nDivSamp);
            for k2 = 1:nDivSamp
                [randNum randIndex] =sort(rand(1,Numbers.samp)); %#ok<NCOMMA,ASGLU>
                oldID(k2,:) = find( (randIndex <= (Numbers.samp+1)*k2/(nDivSamp))&(randIndex > (Numbers.samp+1)*(k2-1)/nDivSamp) );
                X_k_subCL(:,:,k2) = X(:,(randIndex <= (Numbers.samp+1)*k2/nDivSamp) & (randIndex > (Numbers.samp+1)*(k2-1)/nDivSamp),1);
                fitIndex_subCL(k2,:) = fitIndex( (randIndex <= (Numbers.samp+1)*k2/nDivSamp) & (randIndex > (Numbers.samp+1)*(k2-1)/nDivSamp) );
                Y_k_subCL(:,:,k2) = Y_k_CL(1:Numbers_goal, (randIndex <= (Numbers.samp+1)*k2/nDivSamp) & (randIndex > (Numbers.samp+1)*(k2-1)/nDivSamp) );
            end
        else
            inf_indexC = cell(1,1);
            X_k_subCL  = X(:,:,1);
            Y_k_subCL  = Y_k_CL(:,:);
            fitIndex_subCL(1,:) = fitIndex;            
        end

        for k2 = 1:nDivSamp
            valid_idx      = not(isinf(sum(Y_k_subCL(:,:,k2),1)) | isnan(sum(Y_k_subCL(:,:,k2),1)) );
            valid_idx2     = find(valid_idx | fitIndex_subCL(k2,:));
            not_valid_idx2 = find(not(valid_idx) | not(fitIndex_subCL(k2,:)));
            l_not_valid(k2+1) = size(not_valid_idx2,2); 
            l_valid(k2+1)     = size(valid_idx2,2);

            X_k_subCL_vl = X_k_subCL(:,valid_idx2,k2);
            Y_k_subCL_vl = Y_k_subCL(:,valid_idx2,k2);

            X_temp  = [X_k_subCL_vl; ones(1,l_valid(k2+1))];
            A_temp  =  Y_k_subCL_vl-([A(:,:,end,k2),ones(Numbers_goal,1)]*X_temp);   %%
            A_temp(A_temp > 1e+15)  = 1e+15;       %%preventing inf calc.
            A_temp(A_temp < -1e+15) = -1e+15;       %%preventing -inf calc.
            A_temp2 = (pinv(X_temp')*(A_temp)')'...
                 +[A(:,:,end,k2),ones(Numbers_goal,1)];
            A(:,:,end+1,k2) = A_temp2(:,1:Numbers.xfree);

            y_star = y_star_set([1,2,3], 1);
            y_star = y_star.*(1+rand(size(y_star))/0.5*0.1);  % addition of a positive noize aropund target condition (5%)
            
            if nDivSamp>1
                Y_star = repmat(y_star,1,Numbers.samp/nDivSamp );
            else
                Y_star = repmat(y_star,1,Numbers.samp);
            end

            % Replacing notcalc_index2 parameter sets with linear combination of successful samples
            XforCalc = X_k_subCL(:,:,k2);
            if size(inf_indexC{k2},2) ~=0
                X_valid = XforCalc(:,valid_idx2);
                rand_val=rand(l_valid(k2+1),l_not_valid(1));
                rand_val_normalized = rand_val * diag(1./sum(rand_val));
                Xreplace = (X_valid * rand_val_normalized);
                XforCalc(:,not_valid_idx2) = Xreplace;
            end
            if nDivSamp>1
                S(:,:)=A(:,:,end,k2)'* (pinv(((A(:,:,end,k2)*A(:,:,end,k2)'))) *(Y_star-A_temp2*[XforCalc; ones(1,size(XforCalc,2))]));
            else
                S(:,:)=A(:,:,end,k2)'* pinv(((A(:,:,end,k2)*A(:,:,end,k2)'))) *(Y_star-A_temp2*[XforCalc; ones(1,Numbers.samp)]);
            end
            if isnan(A(:,:,1,k2)*A(:,:,1,k2)')
                'error: RCOND equals NaN', break %#ok<NOPRT>
            end 
            X_k_subCL(:,:,k2) = S + XforCalc;
            
            X_mat_tmp = S+XforCalc;
            
            % damping process addition
            k_damp=1;
            while k_damp < 6
                damp_index = (sum( ...
                    (X_mat_tmp > repmat(log(x_uplim),1,size(X_mat_tmp,2))) ...
                    +(X_mat_tmp < repmat(log(x_lowlim),1,size(X_mat_tmp,2))) ...
                    + isnan(X_mat_tmp) + isinf(X_mat_tmp) ) ~=0);
                if sum(damp_index)==0
                    break
                end
                X_mat_tmp(:,damp_index) = X_mat_tmp(:,damp_index)-((1/2)^k_damp).*S(:,damp_index);
                k_damp = k_damp+1;            
            end
            X_k_subCL(:,:,k2) = X_mat_tmp;
        end

        if nDivSamp>1
            Xtemp = [];
            for k2=1:nDivSamp
                Xtemp = [Xtemp X_k_subCL(:,:,k2)];
            end
            X(:,:,1) = Xtemp;
        end

        sumfit = sum(fitIndex);
        if sumfit == 0
            fitsum0_count = fitsum0_count + 1;  % counting a number of no collection in recent loops
        else
            fitsum0_count = 0;
        end
        
        if fitsum0_count > 2 && k_CNM>20
            if k_CNM< 100              % k_CNM < 100, accept 2 sequential loops with no fit
                targetVal = targetVal.*1.03;      %1.03  %target_relaxation when no parameter sets are obtained in recent 2 loops
                fitsum0_count = 0;
            elseif fitsum0_count > 5   % k_CNM >= 100, accept 5 sequential loops with no fit
                targetVal = targetVal.*1.03; %1.03  %target_relaxation when no parameter sets are obtained in recent 5 loops
                fitsum0_count = 0;
                if (sum(duNorm(k_CNM,4)+MRJ(k_CNM,4)+FPJ(k_CNM,4))*0)~=0   % if O_functions return invalid value
                    X(:,:,1) = xmat;
                end
            end
        end
        
        % MAX MIN force limitation
        X_UPLIM  = repmat(x_uplim, 1, size(X(:,:,1),2));  % note normal scale
        X_LOWLIM = repmat(x_lowlim, 1, size(X(:,:,1),2)); % note normal scale
        not_valid_idx2 = find(sum(exp(X(:,:,1)) > X_UPLIM,2) | sum(exp(X(:,:,1)) < X_LOWLIM,2));
        valid_idx2 = find(sum(exp(X(:,:,1)) <= X_UPLIM,2) & sum(exp(X(:,:,1)) >= X_LOWLIM,2)); 

        % for parameter with in range, check current max and min for each set
        % assuming not valid other parameter may be in the same range for
        % each set
        Max_X = repmat(max( X(valid_idx2,:,1) ), size(X,1) ,1);  % max for each set
        Min_X = repmat(min( X(valid_idx2,:,1) ), size(X,1) ,1);  % min for each set        
        Max_X(Max_X > log(X_UPLIM) ) = log(X_UPLIM(Max_X > log(X_UPLIM)));
        Min_X(Min_X < log(X_LOWLIM) ) = log(X_LOWLIM(Min_X < log(X_LOWLIM)));        
        if size(Max_X,1)~=size(X,1)
            Max_X = repmat(log(x_uplim), 1,size(X,2)); 
        end
        if size(Min_X,1)~=size(X,1)
            Min_X = repmat(log(x_lowlim), 1,size(X,2)); 
        end

        X_random_gen_in_range = (Max_X - Min_X).*rand(size(X)) + Min_X;
        X(not_valid_idx2,:,1) = X_random_gen_in_range(not_valid_idx2,:,1);
        
        % Prepare X for NEXT calculation
        XforCalculation = prepareXforCalculation(Numbers,InitialXY,ParameterIndex,X(:,:,1));
        k_CNM = k_CNM +1;
    end
    
    
    %% L-BFGS loop
    % L-bfgs prep
    disp([char(10),char(10),'globalization process ', num2str((k_SS+1)/2),char(10)])

    options_rep = [];
    options_rep.display = 'none';
    options_rep.maxFunEvals = 100;   % 100
    options_rep.maxItr = 100;        %100
    options_rep.optTol = 1e-15;
    options_rep.Method = 'lbfgs';
    options_rep.progTol = 1e-20;
    options_rep.useMex = 0;      
    
    options = [];
    options.display = 'none';
    options.maxFunEvals = 1000;   % default1000;
    options.maxItr = 500;         % default500;
    options.optTol = 1e-15;       % 1e-5;
    options.Method = 'lbfgs';     % 'lbfgs';
    options.progTol = 1e-20;      % 1e-8
    options.useMex = 0;           % mex
    x_ori2 = X(:,:,1);
    tanMaxLength = Numbers.maxTan_shift;  % degree of tangential noise in globalized LBFGS
    ext = Numbers.noise;                  % degree of noise in expanding of obtained sets
    
    parfor_progress(Numbers.samp);
    parfor k1=1:Numbers.samp %#ok<*PFUNK> 
        x = X(:,k1,1);
        lb = log(x_lowlim);
        ub = log(x_uplim);
        option_fmincon = optimoptions(@fmincon,'Display','off');
        warning('off','all')

        objFun2 = @(y) ObjFun(y, Numbers, InitialXY, ParameterIndex, u2, u2ori, uNOTconst, N_obsP, NZLM_target, odeF2, JacbF2,cal_tol);
        objFunNorm2_vec = @(y) ObjFun_dU_vec(y, Numbers, InitialXY, ParameterIndex, u2ori, uNOTconst,odeF2,  cal_tol);
%         objFunNorm2 = @(y) ObjFunNorm(y, Numbers, InitialXY, ParameterIndex, u2ori, uNOTconst, odeF2, cal_tol); %can be used in some cases

        try
            if k_SS==1
                repN = 20;    % noumber of iteration of gL-BFGS (k_SS =1)
            else
                repN = 5;     % noumber of iteration of gL-BFGS (k_SS >1)
            end
            
            % inf or nan correction
            if sum(isnan(x)|isinf(x)) == length(x)
                x = rand(size(x)).*(log(x_uplim)-log(x_lowlim))+ log(x_lowlim);
            end
            x(isinf(x)) = rand(size( x(isinf(x)))).*(log(x_uplim(isinf(x)))-log(x_lowlim(isinf(x))))+ log(x_lowlim(isinf(x)));
            x(isnan(x)) = rand(size( x(isnan(x)))).*(log(x_uplim(isnan(x)))-log(x_lowlim(isnan(x))))+ log(x_lowlim(isnan(x)));
            RetVal = minFuncG_repeat(objFun2, [], x,options_rep,tanMaxLength,repN,ext,x_lowlim,x_uplim,x_ori2);

            outRange = (not(exp(RetVal)<x_uplim)|not(exp(RetVal)>x_lowlim));
            if sum(outRange)~=0
                randReplace = log(x_lowlim) + (log(x_uplim)-log(x_lowlim)).*rand(size(RetVal)); 
                RetVal(outRange) = randReplace(outRange);
            end

            RetVal2 = fsolve(objFunNorm2_vec,RetVal,optimoptions('fsolve','Display','off', 'OptimalityTolerance',1e-6,'StepTolerance',1e-6));
            
            if Numbers.SSS == 0
                RetVal2 = minFunc_global(objFun2,RetVal2,options,0,x_lowlim,x_uplim);      % 3 objs optim. not global, can be used in some cases
            else
                RetVal2 = fmincon(objFun2,RetVal2,[],[],[],[],lb,ub,[],option_fmincon);
            end
            
            X(:,k1,1) = (RetVal2(1:end));
            if sum(isnan(X(:,k1,1)))~=0
                RetVal2 = fmincon(objFun2,x,[],[],[],[],lb,ub,[],option_fmincon);
                X(:,k1,1) = (RetVal2(1:end));
            end
            if sum(isinf(X(:,k1,1)))~=0
                RetVal2 = fmincon(objFun2,x,[],[],[],[],lb,ub,[],option_fmincon);
                X(:,k1,1) = (RetVal2(1:end));
            end
        catch
            RetVal2 = fmincon(objFun2,x,[],[],[],[],lb,ub,[],option_fmincon);
            X(:,k1,1) = (RetVal2(1:end));
        end

        parfor_progress;
    end
    parfor_progress(0);
    Teaps_output_raw(:,:,2*k_SS) = X(:,:,1);
    warning('on','all')
    
    x_store_final = X(:,:,1) ;
    XforCalculation = prepareXforCalculation(Numbers,InitialXY,ParameterIndex,x_store_final);
    save([resultFolder, 'tmpData',num2str(k_SS),'.mat'])
end
delete(gcp)
Teaps_output = InitialXY;
Teaps_output.xfree_distribution = x_store_final;
Teaps_output.X_Teaps = Teaps_output_raw;
Teaps_output.finalX = XforCalculation';
save([resultFolder, 'Teaps_output.mat'],'Teaps_output')


%% Addition of path to optimization module
rmpath(genpath('minFunc_2012'));

end

function [f df] = ObjFun(y, Numbers, InitialXY, ParameterIndex, u2, u2ori, uNOTconst, N_obs, NZLM_target, odeF, JacbF, cal_tol)
    hstep = 1e-6;   % numerical calc of derivertive, the size of slice
    w1 = Numbers.weight_fix;   % 2 (usually fixed)
    w2 = Numbers.weight_relax; % usually values between 2 and 1 will work
    w3 = Numbers.weight_basin; % usually set the value same with w2.
    rj_upper_lim = 1e10;
    X  = prepareXforCalculation2(Numbers,InitialXY,ParameterIndex,y) ;   
    u3ori = u2ori(:,1); %#ok<PFBNS>

    % Calc. O_fix
    dU = odeF(X,u3ori);
    O_fix = norm(dU(uNOTconst),2);
    O_FIX = (O_fix-0)/cal_tol;
%     O_FIX = (O_fix-0);
    
    % Calc. O_relax
    Jori = JacbF(X, u3ori);
    try
        rjori_temp = real(eig(Jori(uNOTconst, uNOTconst)));
    catch
        rjori_temp = rj_upper_lim;
    end
     O_relax = max([(max(rjori_temp(abs(rjori_temp)> cal_tol))-(NZLM_target))/cal_tol,0]);
%      O_relax = exp(max([(max(rjori_temp(abs(rjori_temp)> cal_tol))-(NZLM_target))/cal_tol,0]));
%      O_relax = exp(max([(max(rjori_temp(abs(rjori_temp)> cal_tol))-(NZLM_target))/abs(NZLM_target),0]));
%      O_relax = max([(max(rjori_temp(abs(rjori_temp)> (cal_tol *1e-2) ))-(NZLM_target))/(cal_tol *1e-2),0]);
%      O_relax = max([(max(rjori_temp(abs(rjori_temp)> cal_tol))-(NZLM_target))/abs(NZLM_target),0]);
    if isempty(O_relax)
        O_relax = rj_upper_lim;
    end
    O_RELAX = O_relax;
%     O_RELAX = log(1 + O_relax);  % may be stable in complicated problems
    
    % Calc. O_basin
    O_basin = -Inf;   % initializaiton
    for ks = 1:N_obs
        u3 = u2(:,ks); %#ok<PFBNS>
        J = JacbF(X, u3);
        try
            J_notConst = J(uNOTconst, uNOTconst);
            mrj_tmp = max(real(eig( J_notConst + J_notConst' )));             
        catch
            mrj_tmp = rj_upper_lim;  %tenatative
        end
        O_basin = max(mrj_tmp,O_basin);
        if O_basin > rj_upper_lim
            O_basin = rj_upper_lim;
        end
    end
%     O_BASIN = (O_basin>0).*O_basin;
    O_BASIN = ((O_basin>0).*O_basin)/cal_tol;
%     O_STAB = log(1+(O_basin>0).*O_basin);  % tyr for complicated models

    f = O_FIX^w1 + O_BASIN^w2 + O_RELAX^w3 ;    % obj fun definition

    % calculation of derivertives
    ylength = length(y);
    df=zeros(size(y));
    for k=1:ylength
        y2=y;
        y2(k) = y(k)+hstep;
        X2  = prepareXforCalculation2(Numbers,InitialXY,ParameterIndex,y2);
        u3ori = u2ori(:,1); %#ok<PFBNS>
        
        % Calc. O_fix2
        dU2 = odeF(X2,u3ori);
        O_fix2 = norm(dU2(uNOTconst),2);  %  --> equals to zero L_inf norm
        O_FIX2 = (O_fix2-0)/cal_tol;
%         O_FIX2 = (O_fix2-0);
        
        % Calc. O_relax2
        Jori = JacbF(X2, u3ori);
        try
            rjori_temp = real(eig(Jori(uNOTconst, uNOTconst)));
        catch
            rjori_temp = rj_upper_lim;
        end
        O_relax2 = max([(max(rjori_temp(abs(rjori_temp)> cal_tol))-(NZLM_target))/cal_tol,0]);
%         O_relax2 = exp(max([(max(rjori_temp(abs(rjori_temp)> cal_tol))-(NZLM_target))/cal_tol,0]));
%         O_relax2 = exp(max([(max(rjori_temp(abs(rjori_temp)> cal_tol))-(NZLM_target))/abs(NZLM_target),0]));
%         O_relax2 = max([(max(rjori_temp(abs(rjori_temp)>  (cal_tol *1e-2)))-(NZLM_target))/ (cal_tol *1e-2),0]);
%         O_relax2 = max([(max(rjori_temp(abs(rjori_temp)> cal_tol))-(NZLM_target))/abs(NZLM_target),0]);
        if isempty(O_relax2)
            O_relax2 = rj_upper_lim;
        end
        O_RELAX2 = O_relax2;
%         O_RELAX2 = log(1 + O_relax2);  % may be stable in complicated problems
        
        % Calc. O_basin2
        O_basin2 = -Inf;
        for ks=1:N_obs
            u3 = u2(:,ks); %#ok<PFBNS>
            J = JacbF(X2, u3);
            try
                J_notConst = J(uNOTconst, uNOTconst);
                mrj_tmp = max(real(eig( J_notConst + J_notConst' ))); %--> less than zero                            
            catch
                mrj_tmp = rj_upper_lim; 
            end
            O_basin2 = max(mrj_tmp,O_basin2);
            if O_basin2 > rj_upper_lim
                O_basin2 = rj_upper_lim;
            end
        end
%         O_BASIN2 = (O_basin2>0).*O_basin2;
        O_BASIN2 = ((O_basin2>0).*O_basin2)/cal_tol;
%         O_STAB2 = log(1+(O_stab2>0).*O_stab2);  % for complicated models
        
        f2 = O_FIX2^w1 + O_BASIN2^w2 + O_RELAX2^w3;        
        df(k) = (f2-f)/hstep;   % gradient
    end
end


function [f df] = ObjFunNorm(y, Numbers, InitialXY, ParameterIndex,  u2ori, uNOTconst, odeF, cal_tol)
    hstep = 1e-6;   % numerical calc of derivertive, the size of slice
    X  = prepareXforCalculation2(Numbers,InitialXY,ParameterIndex,y);
    u3ori = u2ori(:,1); %#ok<PFBNS>
    dU = odeF( X,u3ori);
    O_fix = norm(dU(uNOTconst),2);
    O_FIX = (O_fix-0)/cal_tol;
    f = O_FIX;        % obj fun definition
    
    ylength = length(y);
    df=zeros(size(y));
    for k=1:ylength
        y2=y;
        y2(k) = y(k)+hstep;
        X2  = prepareXforCalculation2(Numbers,InitialXY,ParameterIndex,y2);
        u3ori = u2ori(:,1); %#ok<PFBNS>
        dU2 = odeF(X2,u3ori);
        O_fix2 = norm(dU2(uNOTconst),2); 
        O_FIX2 = (O_fix2-0)/cal_tol;
        f2 = O_FIX2 ;        
        df(k) = (f2-f)/hstep;     % gradient
    end
end

% ObjFun_dU_vec(y, Numbers, InitialXY, ParameterIndex, u2ori, uNOTconst, odeF2, cal_tol);
function f = ObjFun_dU_vec(y, Numbers, InitialXY, ParameterIndex, u2ori, uNOTconst, odeF,cal_tol)
% returns dU vector
    X  = prepareXforCalculation2(Numbers,InitialXY,ParameterIndex,y);
    u3ori = u2ori(:,1); %#ok<PFBNS>
    dU = odeF( X,u3ori);
    O_fix = dU(uNOTconst);
    O_FIX = (O_fix-0)./cal_tol;
    f = O_FIX;        % obj fun definition
end


function XforCalculation = prepareXforCalculation(Numbers,InitialXY,ParameterIndex,X_subCL)

XforCalculation = zeros(Numbers.samp,Numbers.xfree+Numbers.xfix);
XforCalculation(1:Numbers.samp,ParameterIndex.fix) = repmat(InitialXY.xfix(1:Numbers.xfix),Numbers.samp,1);
XforCalculation(1:Numbers.samp,ParameterIndex.free) = exp(X_subCL(ParameterIndex.free2,1:Numbers.samp))';
end

function XforCalculation = prepareXforCalculation2(Numbers,InitialXY,ParameterIndex,X_subCL)

XforCalculation = zeros(1,Numbers.xfree+Numbers.xfix);
XforCalculation(1,ParameterIndex.fix) = repmat(InitialXY.xfix(1:Numbers.xfix),1,1);
XforCalculation(1,ParameterIndex.free) = exp(X_subCL(ParameterIndex.free2,1));
end

function u = u_constraint(u, MB_constraint2)
% reflecting mass balance constraints in Model folder
% u: vector, parameter number x 1
%    or matrix, parameter number x number of parameter set

%% applying the constraint
if size(u,2)==1
    u = MB_constraint2(u);
elseif size(u,2)>1
    for k_ps = 1:size(u,2)
        u(:,k_ps) = MB_constraint2(u(:,k_ps));
    end
else
    error('Error in Mass Balance constraint.')
end

end