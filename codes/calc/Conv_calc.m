function [dU_norm NZLM] = Conv_calc( uinitial, xvalue, Model, uconst) %#ok<NCOMMA>
% convergence calc

warning('off','all');

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
cd(curDir)

%%
winitial = uinitial(unoconst);
RHS_ode = @(t, y) ( (odeF2(xvalue,u_ret(uinitial,unoconst,y)))*unoconstM'   );
% RHS_ode = @(t, y) ( (odeF2(xvalue,u_ret(uinitial,unoconst,y)))   );

% norm of dU vector: to see convergence to the fixed point
dU_norm = norm(RHS_ode(0, winitial));

JacFUN = @(t, y, FY) (unoconstM * JacbF2(xvalue,u_ret(uinitial,unoconst,y))* unoconstM' );
jac_fix = JacFUN(0, winitial ,0);
Lambda = real( eig(jac_fix));

% non zero lambda max (NZLM)
NZLM =  max(Lambda(abs(Lambda)> 1e-8 ));

end

function u = u_ret(uinitial,unoconst,w)
    u = uinitial;
    u(unoconst) = w;
end
