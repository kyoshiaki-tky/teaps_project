function sim_fitTC( Teaps_output, Model,u_info, TmaxRun, n_Tdiv , Rth, n_Sim, resultFolder)
%   Based on parameters obtained by the CNM, the kinetic behavior will be
%   simulated

homeDir = pwd;

%% initialization
uconst = u_info.uconst;
uinitial = u_info.uinitial;
uConv_idx = find(u_info.uConv);


%% Calculation
cd('../calc');

fit_idx = Teaps_output.fit_idx;
Xvalue = Teaps_output.finalX(:,fit_idx);

% select  parameter set randomly
nPS = size(Xvalue,2);
if n_Sim >= nPS
    n_Sim = nPS;
    pSetList = 1:n_Sim;
else
    pSetList = randperm(nPS,n_Sim); 
end
        
u_sim_length = length(uConv_idx);
rowN = min([u_sim_length,8]);
colN = ceil(u_sim_length/rowN);

% plot prep
scrsz=get(0, 'ScreenSize');
f1=figure('Position', [1,scrsz(4), scrsz(3)/2.2,scrsz(4)]);
    
for j=1:u_sim_length
    ax1(j)=subplot(rowN,colN,j);
end

% solve ode
for kc2=1:n_Sim
    x_1set = Xvalue(:,pSetList(kc2));  %temp
    uinitial2 = uinitial.*(1+(rand(size(uinitial))-0.5)/0.5*Rth ) ;
    % note that Mass balance constraint is implemented in 'solveCVode'
    
    [t, timeCourse] = solve_CVode(uinitial2, x_1set, Model, uconst, TmaxRun, n_Tdiv);
    for j= 1:u_sim_length
        if kc2>1
            hold on
            set(ax1(j),'NextPlot','add')
        end
        
        plot(ax1(j),t, timeCourse(:,uConv_idx(j)))
        xlabel(ax1(j),'time')
        ylabel(ax1(j),['u' num2str(uConv_idx(j))])   % for constructing article, "u" is substituted by "x"
        xlim(ax1(j),[0 TmaxRun])
        if uconst(uConv_idx(j))
            ylim(ax1(j),[0 1.5*uinitial(uConv_idx(j))])
        end
        if ~uconst(uConv_idx(j))
            ylim(ax1(j),[0 inf])
        end
        yl= ylim(ax1(j));
        if yl(2) > 4*uinitial(uConv_idx(j))
            ylim(ax1(j),[0 1.5*uinitial(uConv_idx(j))])
        end
        if kc2 > 1
            hold off
        end
    end
    hold on
    set(f1,'NextPlot','add')
end

% data save
saveas(f1,[resultFolder,'../simulationOverlay_teaps'],'fig');
saveas(f1,[resultFolder,'../simulationOverlay_teaps'],'epsc');
save([resultFolder,'simulationData.mat'])
cd(homeDir)
end
