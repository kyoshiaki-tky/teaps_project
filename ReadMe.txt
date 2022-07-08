ReadMe for TEAPS (Ver. 1.1 Release)
Only checked in Mac and Matlab 2019 environment.

1. Model Preparation
REQUIREMENT: "CellDesigner", "Maxima"

"CellDesigner" is available from URL below:
http://www.celldesigner.org

"Maxima" is available from homebrew or URL below:
https://maxima.sourceforge.io/
Before execution, the location of "maxima" installed folder need to be specified.

1-1. Prepare a model using CellDesigner XML file.
* Sample models are stored in "SBML_model" folder.
* Need to input kinetic raw for each reaction.
* Currently, description by four arithmetic operations is only acceptable for further process.
* Need to input initial amounts of state variable (species). In this step, just input any number of interset. 
* Need to input kinetic parameter values. In this step, just input any number of interest.
* Saving location does not matter. Save in your favorite folder.

1-2. Converting a model to Matlab format
* If same model is previously analyzed, the model name is renewed.
1-2-1. Open Matlab and change current directory to "tepas_project" folder. 
1-2-2. Execute "modelFilesPrep". This will automatically prepare model folder for tepas.
1-2-3. In the popup dialog, select a model file (CellDesigner XML file).
1-2-3. Input location of directory where "maxima" is installed.
1-2-4. The location of Model Folder and Model Name are displayed for notice. 
1-2-5. Proceed to manual preparation steps shown as further requirements shown in the command windows.
=========================================================
Manually edit 'Parameters.csv' file in model folder.      
   --> Define initial parameter range and search range
Manually edit 'State_variables.csv' file in model folder. 
   --> Define a target state of variables
If required, edit 'MB_constraint.m' file in model folder.  
   --> Define mass balance constraints if required.
     Also see "MB_constraint_example.txt"
Manually edit TEAPS setting in 'main.m'.                  
   --> Able to modify TEAPS setting
     "Model name", "number of parameter sets"
     "target basin size", "a standard of relaxation speed"
      and etc.
=========================================================


2. Execution TEAPS main routine.
REQUIREMENT: "Matlab Parallel Computing toolbox", "Matlab optimization toolbox", "SUNDIALS CVODE solver TB"
SUNDIALS CVODE solver TB" is only available upto version 2.6.0. (In revcent version, matlab toolbox is not included.)
For detail, see the URL below:
https://computing.llnl.gov/projects/sundials

2-1. Execute "main".
2-2. If successfully completed, the number of optimized parameter sets is displayed.
2-3. If successfully completed, time course curves simulated by some optimized parameter sets are displayed.
2-4. Obtained parameter sets are saved in results folder, which is displayed in the end of main routine.
* Some results are shown as examples.


NOTE that codes include parfor_progress.m and minFunc2012, which are made by other people than our team and public in the web.
parfor_progress.m is available from
https://jp.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor-progress-bar-that-works-with-parfor

minFunc2012 is available from
https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html
We modified codes of l-bfgs part for parameter search globalization.

See "LICENSE.txt" for license information. 