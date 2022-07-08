function [t u] = solve_CVode( uinitial, xvalue, Model, uconst, Tmax, Tdiv)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
warning('off','all');

global RHS_ode JacFUN

unoconst = ~uconst;
idx_unoconst = find(unoconst);
unoconstM = zeros(sum(unoconst),length(unoconst));
for row = 1:sum(unoconst)
    unoconstM(row,idx_unoconst(row))=1;
end

%% retrieve ode informaiton
curDir = pwd;
cd([pwd, '/../../Model/', Model])
odeF2= @odeSBML;
JacbF2 = @JacbSBML;
MB_constraint2= @MB_constraint;
cd(curDir)

%% Application of Mass balance constraint against initial state
uinitial = MB_constraint2(uinitial);

%%
TimePoints = 0:(Tmax/(Tdiv-1)):Tmax;
winitial = uinitial(unoconst);

try
    RHS_ode = @(t, y) ( (odeF2(xvalue,u_ret(uinitial,unoconst,y)))*unoconstM'   );
    JacFUN = @(t, y, FY) (unoconstM * JacbF2(xvalue,u_ret(uinitial,unoconst,y))* unoconstM' );
    optionsCV = CVodeSetOptions('LinearSolver','Dense','JacobianFn',@CV_JacFUN, 'MaxNumSteps', 500, 'RelTol',1e-4, 'AbsTol', 1e-6);
    CVodeInit(@CV_RHS_ode,'BDF','Newton', TimePoints(1) , winitial, optionsCV);
    i = 1;
    tt = 0;
    t = 0;   %originally t(i) =0;
    w = winitial;   % originally w(i)=winitial

    while tt < Tmax
        [status,tt,ww] = CVode(Tmax,'OneStep');  %SOLVING ODE
        i = i+1;
        t(:,i) = tt;
        w(:,i) = ww;
    end

    tsize = length(t);
    u = repmat(uinitial', tsize, 1);
    u(:,unoconst) = w';

catch
    t = 100;
    tsize = length(t);
    u = repmat(uinitial', tsize, 1);
    u(:,unoconst) = 100;  % to lead not SBS
end


end

%%

function u = u_ret(uinitial,unoconst,w)
    u = uinitial;
    u(unoconst) = w;
end

function [hReturn flag]=CV_RHS_ode(t,w)
global RHS_ode
flag = zeros(size(w));
hReturn = RHS_ode(t,w);

end
function [J, FLAG] = CV_JacFUN(t, w, FY)
%The input argument FY contains the current valu of f(t,y)
%  In this case FY corresponds to hret maybe
global JacFUN

FLAG=zeros(size(w));
J = JacFUN(t, w, FY);
end
%