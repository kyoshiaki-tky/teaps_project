function [uconst uinitial uConv] = importFuncSBML( Model , varargin) %#ok<NCOMMA>
%Extraction of ODE flux from SBML

Address = ['../../Model/',Model, '/', Model, '.xml'];
fileID = fopen(Address);
SBMLc = textscan(fileID,'%s','delimiter','\n');
fclose(fileID);

cur = pwd;
cd(['../../Model/',Model, '/'])
if exist([Model, '_Ans_SS.txt'],'file')||~isempty(varargin)
    Address2 = ['../../Model/',Model, '/', Model, '_Ans_SS.txt'];
    fileID2 = fopen(Address2);
    a2x = textscan(fileID2,'%s','delimiter','=,');
    fclose(fileID2);

    a2x = a2x{1}(~cellfun(@isempty, a2x{1}));
    a2x=strrep(a2x,'[','');
    a2x=strrep(a2x,']];','');
    a2x = reshape(a2x,2,[]);
    a2x = a2x';
    a2x=strrep(a2x,' ','');
end
cd(cur)

SBMLc=strrep(SBMLc{1},'celldesigner:class', 'CLASS_type');
SBMLc=SBMLc(cellfun(@isempty, cellfun(@(Y) strfind(Y, 'celldesigner'), SBMLc,'UniformOutput', false)));


RowN = 1;
un=0;
an=0;
u2s={};
globalp=true;

%% sp2u conversion with judgement of conc/amount and global parameters aquisition(x)
while RowN <= length(SBMLc)
    if strfind(SBMLc{RowN},'<species ') %metaid="') 
        icstart = strfind(SBMLc{RowN},'initialConcentration="')+22;
        iastart = strfind(SBMLc{RowN},'initialAmount="')+15;
        if ~isempty(icstart)
            un = un+1;
            dcIndex=[];
            dcIndex = strfind(SBMLc{RowN},'"');
            icend = dcIndex(find(dcIndex==(icstart-1))+1)-1;
            sidstart = strfind(SBMLc{RowN},' id="')+5;
            sidend = dcIndex(find(dcIndex==(sidstart-1))+1)-1;
            namestart = strfind(SBMLc{RowN},' name="')+7;    %added
            nameend = dcIndex(find(dcIndex==(namestart-1))+1)-1;    %added
            notcnst = isempty(strfind(SBMLc{RowN},'constant="true"'));  %added
            u2s(un,:) = {['u(' num2str(un) ')'],SBMLc{RowN}(sidstart:sidend),notcnst, SBMLc{RowN}(icstart:icend), 'concentration',SBMLc{RowN}(namestart:nameend),[]};  %modified
            SBMLc = strrep(SBMLc, u2s(un,2),u2s(un,1));
            dU{un}='0';
        elseif  ~isempty(iastart)
            un = un+1;
            dcIndex = [];
            dcIndex = strfind(SBMLc{RowN},'"');
            iaend = dcIndex(find(dcIndex==(iastart-1))+1)-1;
            sidstart = strfind(SBMLc{RowN},' id="')+5;
            sidend = dcIndex(find(dcIndex==(sidstart-1))+1)-1;
            namestart = strfind(SBMLc{RowN},' name="')+7;    %added
            nameend = dcIndex(find(dcIndex==(namestart-1))+1)-1;    %added
            notcnst = isempty(strfind(SBMLc{RowN},'constant="true"'));  %added
            u2s(un,:) = {['u(' num2str(un) ')'],SBMLc{RowN}(sidstart:sidend),notcnst, SBMLc{RowN}(iastart:iaend), 'amount',SBMLc{RowN}(namestart:nameend),[]};  %modified
            SBMLc = strrep(SBMLc, u2s(un,2),u2s(un,1));
            dU{un}='0';
        end
        RowN = RowN+1;
    elseif strfind(SBMLc{RowN},'CLASS_type')  %steady state judgement  unless class is DEGRADED or PHENOTYPE, the values should be converged
        if un == 0
        elseif ~isempty(strfind(SBMLc{RowN},'DEGRADED'))||~isempty(strfind(SBMLc{RowN},'PHENOTYPE'))
            u2s(un,7) = {'divergent'};  %% may be divergent
            uConv(un) = false;
        else
            u2s(un,7)= {'converged'};  %% should be converged
            uConv(un) = true;
        end
        RowN = RowN+1;
    elseif strfind(SBMLc{RowN},'<compartment ') %metaid="')
        dcIndex3 = strfind(SBMLc{RowN},'"');
        cidstart = strfind(SBMLc{RowN},' id="')+5;
        cidend = dcIndex3(find(dcIndex3==(cidstart-1))+1)-1;
        csizestart = strfind(SBMLc{RowN},'size="')+6;
        csizeend = dcIndex3(find(dcIndex3==(csizestart-1))+1)-1;
        an = an+1;
        a2p(an,:) = {['a(' num2str(an) ')'],SBMLc{RowN}(cidstart:cidend),SBMLc{RowN}(csizestart:csizeend), 'Volume'};
        SBMLc = strrep(SBMLc, a2p(an,2),a2p(an,1));
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
            an = an +1;
            a2p(an,:) = {['a(' num2str(an) ')'],SBMLc{RowN}(ParamNameStart:ParamNameEnd),SBMLc{RowN}(ParamValStart:ParamValEnd), 'global'};
            SBMLc = strrep(SBMLc, a2p(an,2),a2p(an,1));
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
    xnstart=an+1;
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
        if strfind(SBMLc{RowN},'<parameter ') %metaid="')
            dcIndex3 = strfind(SBMLc{RowN},'"');
            ParamNameStart = strfind(SBMLc{RowN},' id="')+5;
            ParamNameEnd = dcIndex3(find(dcIndex3==(ParamNameStart-1))+1)-1;
            ParamValStart = strfind(SBMLc{RowN},'value="')+7;
            ParamValEnd = dcIndex3(find(dcIndex3==(ParamValStart-1))+1)-1;

            an = an +1;
            a2p(an,:) = {['a(' num2str(an) ')'],SBMLc{RowN}(ParamNameStart:ParamNameEnd),SBMLc{RowN}(ParamValStart:ParamValEnd), SBMLc{RxnStartIndex(k)}(ReStart:ReEnd)};
        end
        localRe = [];
        localRe ={SBMLc(lrstart:lrend)};
        localRe{1};
        flux{k}=MML2MatEq2(localRe{1});
    end
    a2ptmp=a2p;
    for j=xnstart:an
        flux(k) = strrep(flux{k}, a2ptmp(j,2),a2ptmp(j,1));
        a2ptmp(xnstart:an,2)=strrep(a2ptmp(xnstart:an,2),a2ptmp(j,2),a2ptmp(j,1));
    end
end

dUSBML=dU;
for k =1:size(flux,2)
    dUSBML= strrep(dUSBML, ['flux(' int2str(k) ')'], flux(k));
end
dUSBML= strrep(dUSBML, ' ','');

uConvList = logical(cellfun(@isempty,cellfun(@(Y) strfind(Y, 'divergent'),u2s(:,7),'UniformOutput', false)));

uconst = ~logical(cell2mat(u2s(:,3)));
uconst = (uconst | not(uConvList));

uinitial = {u2s(:,4)};

end

