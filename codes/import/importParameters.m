function [Numbers InitialXY ParameterIndex u_info] = importParameters(Model, Numbers) %#ok<NCOMMA>
% function [Numbers InitialXY ParameterIndex uinitial uconst] = importParameters(Model, Numbers) %#ok<NCOMMA>
% input flux import etc.
% Initial parameter distributution and Search space are loaded from 'Parameters.csv'
% A part of model information is loaded from 'MODEL.xml' (CellDesigner file)
% Ttarget steady state is loaded from 'State_variables.csv'

%% Common values
InitialXY = struct;  ParameterIndex = struct;

%% Parameter space values

para_source = importStr(['../../Model/',Model, '/','Parameters.csv'], ',');
parameters = str2double((para_source(2:end,4:end))');
ParameterIndex.free = find(parameters(3,:)~=0); Numbers.xfree = size(ParameterIndex.free,2);
ParameterIndex.fix = find(parameters(3,:)==0); Numbers.xfix = size(ParameterIndex.fix,2);
ParameterIndex.free2 = 1:Numbers.xfree ;

Numbers.input = Numbers.xfree + Numbers.xfix; 

% import initial min&max parameters, logarithmic search
InitialXY.xfree_import = log(parameters(1:2,ParameterIndex.free)); 

% import Search space min & max parameters (later: logarithmic conversion)
InitialXY.xfree_ss_import = (parameters(4:5,ParameterIndex.free));

% import fixed parameter value
InitialXY.xfix = parameters(1,ParameterIndex.fix);

%% Importing SBML model
[uconst uinitial_sbml uConv] = importFuncSBML(Model); %#ok<NCOMMA>

%% Importing Target steady state, 
state_source = importStr(['../../Model/',Model, '/','State_variables.csv'], ',');
uinitial = str2double((state_source(2:end,5)));
if prod(size(uinitial)==size(uinitial_sbml{1}))~=1
    error('Error in importing a target steady state.')
end

u_info = struct;
u_info.uinitial = uinitial;
u_info.uconst = uconst;
u_info.uConv = uConv;

end