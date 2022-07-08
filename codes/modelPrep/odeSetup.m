function odeSetup( SBML, modelFolder)
% Preparing 'JacbSBML.m' file' and 'odeSBML.m' file from "maxima" output

Address1 = [modelFolder '/', SBML, '_jacobian.txt'];
fileID1 = fopen(Address1);
jac = textscan(fileID1, '%s');
fclose(fileID1);

Address2 = [modelFolder '/', SBML, '_dUSBML_xu.csv'];
fileID2 = fopen(Address2);
ode = textscan(fileID2, '%s');
fclose(fileID2);


jac = strrep(jac{1}, 'matrix(', 'dUJ=[');
jac = strrep(jac, ');', '];');
jac = strrep(jac, '],','];');
jac2(1,:) = jac(:);
clear jac;
dUJac = strrep(strjoin(jac2), ' ', '');
clear jac2;

ode=ode{1};
usize = length(ode);
for k = 1: usize
    lhs{k,1} =['dUSBML(' int2str(k) ')='];
end
dUSBML = strcat(lhs, ode, ';');

Header1 = ['function dUJ = JacbSBML(x, u)'];
Footer = ['end'];

Header2 = ['function dUSBML = odeSBML(x, u)'];
uNeverNegative = ['u(u<0)=0;'];


filename = [modelFolder, '/JacbSBML.m'];
fid = fopen(filename, 'w');
fprintf(fid, '%s\n', Header1);
fprintf(fid, '%s\n', uNeverNegative);
fprintf(fid, '%s\n', dUJac);
fprintf(fid, '%s\n', Footer);
fclose(fid);

filename2 = [modelFolder, '/odeSBML.m'];
fid2 = fopen(filename2, 'w');
fprintf(fid2, '%s\n', Header2);
fprintf(fid2, '%s\n', uNeverNegative);
for row = 1:usize
    fprintf(fid2, '%s\n', dUSBML{row,1});
end
fprintf(fid2, '%s\n', Footer);
fclose(fid2);

end