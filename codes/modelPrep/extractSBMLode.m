function extractSBMLode( Model, modelFolder, maxima_path)
%Extraction of ODE flux from SBML

Address = [modelFolder, '/', Model, '.xml'];
fileID = fopen(Address);
SBMLc = textscan(fileID,'%s','delimiter','\n');
fclose(fileID);

%Importing simultaneous equation results for steady state.
SBMLc=strrep(SBMLc{1},'celldesigner:class', 'CLASS_type');
SBMLc=SBMLc(cellfun(@isempty, cellfun(@(Y) strfind(Y, 'celldesigner'), SBMLc,'UniformOutput', false)));
RowN = 1;
un=0;
xn=0;
u2s={};
globalp=true;

%% spieces 2 u (state variables in matlab) conversion with judgement of conc/amount and global parameters aquisition(x)
while RowN <= length(SBMLc)
    if strfind(SBMLc{RowN},'<species ') %metaid="') 
        icstart = strfind(SBMLc{RowN},'initialConcentration="')+22;   % concentration based description
        iastart = strfind(SBMLc{RowN},'initialAmount="')+15;          % amount based description
        if ~isempty(icstart)
            un = un+1;
            dcIndex=[];
            dcIndex = strfind(SBMLc{RowN},'"');
            icend = dcIndex(find(dcIndex==(icstart-1))+1)-1;
            sids_tart = strfind(SBMLc{RowN},' id="')+5;
            sid_end = dcIndex(find(dcIndex==(sids_tart-1))+1)-1;
            name_start = strfind(SBMLc{RowN},' name="')+7;    %added
            name_end = dcIndex(find(dcIndex==(name_start-1))+1)-1;    %added
            comp_start = strfind(SBMLc{RowN},' compartment="')+14;    %added
            comp_end = dcIndex(find(dcIndex==(comp_start-1))+1)-1;    %added
            notcnst = isempty(strfind(SBMLc{RowN},'constant="true"'));  %added
            u2s(un,:) = {['u(' num2str(un) ')'],SBMLc{RowN}(sids_tart:sid_end),...
                notcnst, SBMLc{RowN}(icstart:icend), 'concentration',...
                SBMLc{RowN}(name_start:name_end),SBMLc{RowN}(comp_start:comp_end), ...
                []};  %modified
            SBMLc = strrep(SBMLc, u2s(un,2),u2s(un,1));
            dU{un}='0';  % flux initialization
        elseif  ~isempty(iastart)
            un = un+1;
            dcIndex = [];
            dcIndex = strfind(SBMLc{RowN},'"');
            iaend = dcIndex(find(dcIndex==(iastart-1))+1)-1;
            sids_tart = strfind(SBMLc{RowN},' id="')+5;
            sid_end = dcIndex(find(dcIndex==(sids_tart-1))+1)-1;
            name_start = strfind(SBMLc{RowN},' name="')+7;    %added
            name_end = dcIndex(find(dcIndex==(name_start-1))+1)-1;    %added
            comp_start = strfind(SBMLc{RowN},' compartment="')+14;    %added
            comp_end = dcIndex(find(dcIndex==(comp_start-1))+1)-1;    %added
            notcnst = isempty(strfind(SBMLc{RowN},'constant="true"'));  %added
            u2s(un,:) = {['u(' num2str(un) ')'],SBMLc{RowN}(sids_tart:sid_end), ...
                notcnst, SBMLc{RowN}(iastart:iaend), 'amount', ...
                SBMLc{RowN}(name_start:name_end), SBMLc{RowN}(comp_start:comp_end), ...
                []};  %modified
            SBMLc = strrep(SBMLc, u2s(un,2),u2s(un,1));
            dU{un}='0';  % flux initialization
        end
        RowN = RowN+1;
    elseif strfind(SBMLc{RowN},'CLASS_type')  %steady state judgement  unless class is DEGRADED or PHENOTYPE, the values should be converged
        if un == 0
        elseif ~isempty(strfind(SBMLc{RowN},'DEGRADED'))||~isempty(strfind(SBMLc{RowN},'PHENOTYPE'))
            u2s(un,8) = {'divergent'};  %% may be divergent
            uConv(un) = false;
        else
            u2s(un,8)= {'converged'};  %% should be converged
            uConv(un) = true;
        end
        RowN = RowN+1;
    elseif strfind(SBMLc{RowN},'<compartment ') %metaid="')
        dcIndex3 = strfind(SBMLc{RowN},'"');
        cidstart = strfind(SBMLc{RowN},' id="')+5;
        cidend = dcIndex3(find(dcIndex3==(cidstart-1))+1)-1;
        csizestart = strfind(SBMLc{RowN},'size="')+6;
        csizeend = dcIndex3(find(dcIndex3==(csizestart-1))+1)-1;
        xn = xn+1;
        x2p(xn,:) = {['x(' num2str(xn) ')'],SBMLc{RowN}(cidstart:cidend),SBMLc{RowN}(csizestart:csizeend), 'Volume'};            
        SBMLc = strrep(SBMLc, x2p(xn,2),x2p(xn,1));
        RowN = RowN+1;
    elseif strfind(SBMLc{RowN},'<listOfReactions>')
        globalp = false;
        RowN = RowN+1;
    elseif strfind(SBMLc{RowN},'</listOfReactions>')
        globalp = true;
        RowN = RowN+1;
    elseif globalp  %%local parameter rejection
        if strfind(SBMLc{RowN},'<parameter ') %metaid="')
            dcIndex3 = strfind(SBMLc{RowN},'"');
            ParamNameStart = strfind(SBMLc{RowN},' id="')+5;
            ParamNameEnd = dcIndex3(find(dcIndex3==(ParamNameStart-1))+1)-1;
            ParamValStart = strfind(SBMLc{RowN},'value="')+7;
            ParamValEnd = dcIndex3(find(dcIndex3==(ParamValStart-1))+1)-1;
            xn = xn +1;
            x2p(xn,:) = {['x(' num2str(xn) ')'],SBMLc{RowN}(ParamNameStart:ParamNameEnd),SBMLc{RowN}(ParamValStart:ParamValEnd), 'global'};            
            SBMLc = strrep(SBMLc, x2p(xn,2),x2p(xn,1));
        end
        RowN = RowN+1;
    else
        RowN = RowN+1;
    end

end


for k=1:size(u2s,1)
    u2s(:,2) = strrep(u2s(:,2),u2s(k,1),u2s(k,2));
end

%% Local reaction
RxnStartIndex = find(cellfun(@isempty, cellfun(@(Y) strfind(Y, '<reaction'), SBMLc,'UniformOutput', false))==0);
RxnEndIndex = find(cellfun(@isempty, cellfun(@(Y) strfind(Y, '</reaction>'), SBMLc,'UniformOutput', false))==0);
RxnNum = length(RxnStartIndex);
for k = 1: RxnNum
    RowN = RxnStartIndex(k);
    dcIndex = strfind(SBMLc{RowN},'"');
    ReStart = dcIndex(3)+1;
    ReEnd = dcIndex(4)-1;
    flux2re(k,:) = {['flux(' num2str(k) ')'], SBMLc{RxnStartIndex(k)}(ReStart:ReEnd)};
    xnstart=xn+1;
    while RowN <= RxnEndIndex(k)
        if strfind(SBMLc{RowN},'<listOfReactants>')
            RowN = RowN +1;
            while isempty(strfind(SBMLc{RowN},'</listOfReactants>'))
                if strfind(SBMLc{RowN},'<speciesReference')
                    dcIndex2 = strfind(SBMLc{RowN},'"');
                    ReactantStart = dcIndex2(3)+3;
                    ReactantEnd = dcIndex2(4)-2;
                    Unumber = str2num(SBMLc{RowN}(ReactantStart:ReactantEnd));  %added
                    if u2s{Unumber,3}   %added
                        dU{Unumber}=[dU{Unumber}, '-flux(', num2str(k), ')'];
                    end
                end
                RowN = RowN+1;
            end
        elseif strfind(SBMLc{RowN},'<listOfProducts>')
            RowN = RowN +1;
            while isempty(strfind(SBMLc{RowN},'</listOfProducts>'))
                if strfind(SBMLc{RowN},'<speciesReference')
                    dcIndex2 = strfind(SBMLc{RowN},'"');
                    ReactantStart = dcIndex2(3)+3;
                    ReactantEnd = dcIndex2(4)-2;
                    Unumber = str2num(SBMLc{RowN}(ReactantStart:ReactantEnd));  %added
                    if u2s{Unumber,3}   %added
                        dU{Unumber}=[dU{Unumber}, '+flux(', num2str(k), ')'];
                    end
                end
                RowN = RowN+1;
            end
        else
           RowN = RowN+1;
        end
        lrstart = RxnStartIndex(k);
        lrend = RxnEndIndex(k);
        if strfind(SBMLc{RowN},'<parameter metaid="')
            dcIndex3 = strfind(SBMLc{RowN},'"');
            ParamNameStart = strfind(SBMLc{RowN},' id="')+5;
            ParamNameEnd = dcIndex3(find(dcIndex3==(ParamNameStart-1))+1)-1;
            ParamValStart = strfind(SBMLc{RowN},'value="')+7;
            ParamValEnd = dcIndex3(find(dcIndex3==(ParamValStart-1))+1)-1;
            xn = xn +1;
            x2p(xn,:) = {['x(' num2str(xn) ')'],SBMLc{RowN}(ParamNameStart:ParamNameEnd),SBMLc{RowN}(ParamValStart:ParamValEnd), SBMLc{RxnStartIndex(k)}(ReStart:ReEnd)};
        end
        localRe = [];
        localRe ={SBMLc(lrstart:lrend)};
        flux{k} = MML2MatEq(localRe{1});
    end
    X2ptmp = x2p;
    for j=xnstart:xn
        flux(k) = strrep(flux{k}, X2ptmp(j,2),X2ptmp(j,1));
        X2ptmp(xnstart:xn,2)=strrep(X2ptmp(xnstart:xn,2),X2ptmp(j,2),X2ptmp(j,1));
    end
end

x2p_c = size(x2p,2);
x2p2 = x2p;

% set default values
x2p2(:,x2p_c+1) = {0.001};  % initial lower limit
x2p2(:,x2p_c+2) = {1000};   % initial upper limit
x2p2(:,x2p_c+3) = {1};      % 0:fix or 1:free
x2p2(:,x2p_c+4) = {5e-3};      % search area lower limit
x2p2(:,x2p_c+5) = {5e+4};      % search area upper limit
col_rowN = find(strcmp(x2p2(:,4),'Volume'));  % search volume index
x2p2(col_rowN,x2p_c+3) = {0};   % 'Volume' parameter is treated as fixed
for k_v = 1:length(col_rowN)
    x2p2{k_v,x2p_c+1} = cellfun(@str2double,x2p(k_v,3));   % volume is used original value
    x2p2{k_v,x2p_c+2} = cellfun(@str2double,x2p(k_v,3));   % volume is used original value
    x2p2{k_v,x2p_c+4} = cellfun(@str2double,x2p(k_v,3));   % volume is used original value
    x2p2{k_v,x2p_c+5} = cellfun(@str2double,x2p(k_v,3));   % volume is used original value
end

x2p2(:,3)=[];   % omit parameter value described in Cell designer model

x2p2_label={'name in matlab','name in CellDesigner','reaction name',...
    'initial lower limit','initial upper limit','0:fix or 1:free', ...
    'search area lower limit','search area upper limit'};
x2p2 = [x2p2_label; x2p2];


for k_u = 1:size(u2s,1)
    if cell2mat(u2s(k_u,3))
        u2s{k_u,3} = ['not_cnst'];
    else
        u2s{k_u,3} = ['constant'];
    end
end

uConvList = logical(cellfun(@isempty,cellfun(@(Y) strfind(Y, 'divergent'),u2s(:,7),'UniformOutput', false)));
dU(not(uConvList)) = {'0'};
dUSBML=dU;

for k =1:size(flux,2)
    dUSBML= strrep(dUSBML, ['flux(' int2str(k) ')'], flux(k));
end
dUSBML= strrep(dUSBML, ' ','');

% Script for Jacobian calculation with Maxima
odeSet = ['ode:[', strjoin(dUSBML,','), ']$'];
env = ['set_display(''none)$'];
varSet = ['var:[', strjoin(u2s(:,1)',','), ']$'];
output =['jacobian(ode,var)$'];
output2=['stringout("', modelFolder '/',  Model, '_jacobian.txt", %)$'];

[nrows,ncols]= size(dUSBML);
filename = [modelFolder,'/',Model '_dUSBML_xu.csv'];
fid = fopen(filename, 'w');
for row=1:nrows
    fprintf(fid, '%s\n', dUSBML{row,:});
end
fclose(fid);

u2s2 = u2s;
u2s2(:,6)=[];
u2s2_label = {'name in matlab','name in CellDesigner','constant?',...
    'target value','amount/conc','compartment','should be'};
u2s2 = [u2s2_label;u2s2];
u2s2 = u2s2(:,[1,2,6,3,4,5,7]);  % sorting

for k_u = 2:size(u2s2,1)
    x2p_row = (strcmp(u2s2(k_u,3),x2p(:,1)));
    u2s2(k_u,3) = {[u2s2{k_u,3},' : ', x2p{x2p_row,2}]};
end

[nrows,ncols]= size(u2s2);
filename = [modelFolder,'/State_variables.csv'];
fid = fopen(filename, 'w');
fprintf(fid, '%s, %s, %s, %s, %s, %s, %s\n', u2s2{1,:});
for row=2:nrows
    fprintf(fid, '%s, %s, %s, %s, %s, %s, %s\n', u2s2{row,:});
end
fclose(fid);


[nrows,ncols]= size(x2p2);
filename = [modelFolder,'/Parameters.csv'];
fid = fopen(filename, 'w');
fprintf(fid, '%s, %s, %s, %s,%s, %s, %s, %s\n', x2p2{1,:});
for row=2:nrows
    fprintf(fid, '%s, %s, %s, %f,%f, %1.0f, %f, %f\n', x2p2{row,:});
end
fclose(fid);

filename =  [modelFolder,'/',Model '_for_JacCalc.mac'];
fid = fopen(filename, 'w');
fprintf(fid, '%s\n', odeSet);
fprintf(fid, '%s\n', varSet);
fprintf(fid, '%s\n', env);
fprintf(fid, '%s\n', output);
fprintf(fid, '%s\n', output2);
fclose(fid);

system([maxima_path, '/maxima -b ', modelFolder,'/',Model,'_for_JacCalc.mac']);
odeSetup(Model, modelFolder);  % convert jacobian and ode to matlab function

end