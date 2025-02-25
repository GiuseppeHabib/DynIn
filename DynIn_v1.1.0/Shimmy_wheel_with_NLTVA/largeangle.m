function [position,isterminal,direction] = largeangle(t,xt)

value=max([abs(xt(1)),abs(xt(3))])-pi/2;
if value>0
    position=0;
else
    position=1;
end
if max([abs(xt(4)),abs(xt(5))])>100
    position=0;
end
% position=xt(1); % this should mean that we stop when x1=0 (check also this)
isterminal=1;
direction=-1;

end