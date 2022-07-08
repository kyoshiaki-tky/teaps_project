function modelFilesPrep
% REQUIREMENT: maxima installation (and speficy installed molder name)
% select a SBML file to import according to dialog
% input folder where 'maxima' exists.

%% model selection
disp('Select a CellDesigner XML file...')
[model_file,model_path] = uigetfile('*.xml','Select a CellDesigner XML file...');
if isequal(model_file,0)
    error('User selected Cancel')
else
    disp(['User selected ', fullfile(model_path, model_file)])
end

make_modelFolder = [pwd, '/../Model/', model_file(1:end-4)];

%% maxima prep
prompt = ['Which folder maxima exists? :    '];
maxima_path = input(prompt,'s');
if isempty(maxima_path)
    error('specify maxima folder.')
end


%% model folder prep
k = 1;
if ~exist(make_modelFolder, 'dir')
    mkdir(make_modelFolder)
    SBML_name = model_file(1:end-4);
    disp(['Model name: ', SBML_name])    
else
    while exist(make_modelFolder, 'dir')
        k = k+1;
        make_modelFolder = [pwd, '/../Model/', model_file(1:end-4), '_',num2str(k)];
        new_model_file = [model_file(1:end-4), '_',num2str(k),'.xml'];
    end
    mkdir(make_modelFolder)
    disp(['Same name model exists! Model imported as new name.']);
    SBML_name = new_model_file(1:end-4);
    disp(['New Model name: ', SBML_name ]);
end


copyfile([model_path, model_file], [make_modelFolder, '/', SBML_name,'.xml']);
cd('./modelPrep')
copyfile([pwd '/MB_constraint_template.txt'], [make_modelFolder, '/MB_constraint.m']);

%% model ode and jacobian prep 
disp(char(10))
disp('SBML model extraction and Jacobian calculaiton by maxima ...')

extractSBMLode(SBML_name, make_modelFolder,maxima_path);

disp(char(10))
disp('Preparation completed!')
disp(char(10))
disp(['Prepared Model name:  ', char(10), SBML_name ]);
disp(['Prepared Model folder:  ', char(10), make_modelFolder ]);
disp(char(10))
disp('Further requirements:')    
disp('Manually edit ''Parameters.csv'' file in model folder.')    
disp('Manually edit ''State_variables.csv'' file in model folder.')    
disp('If required, edit ''MB_constraint.m'' file in model folder.')    
disp('Manually edit TEAPS setting in ''main.m''. ')    

cd('../')
end