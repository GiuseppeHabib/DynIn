% double well (unforced)
function xprime = sistema(t,x,par)
%par = [zeta a b];
    xprime=zeros(2,1);
    % ============================= State Space ========================
    xprime(1)=x(2);
    xprime(2)=-2*par(1)*x(2)+par(2)*x(1)-par(3)*x(1)^3;
    % ===================================================================
end
