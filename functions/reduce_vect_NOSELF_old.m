function [red_vect,eliminated]=reduce_vect_NOSELF(v1,v2,R,weight)

% dim=length(v(1,:));
n1=length(v1(:,1));
n2=length(v2(:,1));
eliminated=zeros(n1,1);
Rq=R^2;

for i1=1:n1
    for i2=1:n2
        if normweight_q(v1(i1,:)-v2(i2,:),weight)<Rq
            eliminated(i1)=1;
            break
        end
    end
end
red_vect=v1(eliminated==0,:);