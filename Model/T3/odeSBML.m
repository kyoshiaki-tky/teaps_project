function dUSBML = odeSBML(x, u)
u(u<0)=0;
dUSBML(1)=0-(x(2)*u(1))+((x(3)*u(2)*u(3))/(x(4)+u(2)))+x(5);
dUSBML(2)=0+(x(2)*u(1))-((x(3)*u(2)*u(3))/(x(4)+u(2)))-(x(6)*u(2));
dUSBML(3)=0;
dUSBML(4)=0-x(5)+(x(6)*u(2));
end
