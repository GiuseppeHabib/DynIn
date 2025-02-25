function xprime = sistema(t,x,par)
% chain of 4 masses (in modal coordinates)
% modal coordinates are obtained analytically. The derivation of the modal
% coordinates is not discussed here. For details about the mechanical
% system look at % Habib, G. (2021). Dynamical integrity assessment of 
% stable equilibria: a new rapid iterative procedure. Nonlinear Dynamics, 
% 106(3), 2073-2096.

xprime=zeros(8,1);
% ============================= State Space ========================
x2=par.vect(2,1)*x(1)+par.vect(2,2)*x(3)+par.vect(2,3)*x(5)+par.vect(2,4)*x(7);
b=(-(1/6))*sqrt(-289 + 144*par.lr^2) + (4913*x2)/(864*par.lr^2) + (2*par.lr*(sqrt(-289 + 144*par.lr^2) - 17*x2))/sqrt(144*par.lr^2 + 17*x2*(-2*sqrt(-289 + 144*par.lr^2) + 17*x2));
damp=par.Cs*[x(2); x(4); x(6); x(8)];
xprime(1)=x(2);
xprime(2)=-damp(1)-par.val(1,1)*x(1)-par.Ut(1,2)*b;
xprime(3)=x(4);
xprime(4)=-damp(2)-par.val(2,2)*x(3)-par.Ut(2,2)*b;
xprime(5)=x(6);
xprime(6)=-damp(3)-par.val(3,3)*x(5)-par.Ut(3,2)*b;
xprime(7)=x(8);
xprime(8)=-damp(4)-par.val(4,4)*x(7)-par.Ut(4,2)*b;
% ===================================================================
end
