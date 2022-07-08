function matEq = MML2MatEq(MML)
% MML is Eq written in MathML, sequential cells including strings
% Bug fix for reaction starting with "-" @181129

matEq = '';
RowN = 1;
opN = 0;
while RowN <= length(MML)
    if strfind(MML{RowN},'<apply>')
        matEq = [matEq '('];
        RowN = RowN+1;
        opN = opN+1;
        MML{RowN};
        OP(opN) = operatorJudge(MML{RowN});
        countOP_use(opN) = 0;    %%
        RowN = RowN+1;
    elseif strfind(MML{RowN},'<c')
        IntMed = fliplr(strtok(fliplr(strtok(strtok(MML{RowN},'<'),'<')),'>'));
        matEq = [matEq IntMed];
        RowN = RowN+1;
        if strfind(MML{RowN}, '</math>')
            RowN=RowN+1;
        elseif isempty(strfind(MML{RowN},'</apply>'))
            matEq = [matEq OP(opN)];
            countOP_use(opN)=countOP_use(opN)+1;    %%
        end
    elseif strfind(MML{RowN},'</apply>')
        if countOP_use(opN)==0                                 %%
            if strcmp(OP(opN),'-') && (strcmp(OP(opN-1),'*')||strcmp(OP(opN-1),'/'))    %%
                matEq = [matEq '*(-1)'];                       %%
            end                                                %%
        end                                                    %%
        opN = opN-1;
        matEq = [matEq ')'];
        RowN = RowN+1;
        if strfind(MML{RowN},'<apply>')
            matEq = [matEq OP(opN)];
            countOP_use(opN)=countOP_use(opN)+1;    %%
        elseif strfind(MML{RowN},'<c')
            matEq = [matEq OP(opN)];
            countOP_use(opN)=countOP_use(opN)+1;    %%
        elseif strfind(MML{RowN},'</math>')
            RowN=RowN+1;
        end
    else
        RowN = RowN+1;
    end
end 


end
function OC = operatorJudge(X)
if strfind(X,'<plus/>')
    OC = '+';
elseif strfind(X,'<minus/>')
    OC = '-';
elseif strfind(X,'<times/>')
    OC = '*';
elseif strfind(X,'<divide/>')
    OC = '/';
elseif strfind(X,'<power/>')
    OC = '^';
end
end
