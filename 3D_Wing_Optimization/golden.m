%% Golden Section Search


function [alpha]=golden(x, gradf,n)

%% Golden Section Search Parameters to det. min alpha
tau=0.381966;

lb=0;   % Lower Bound Alpha
ub=1;   % Upper Bound Alpha

s=-gradf;

alpha1=(1-tau)*lb+tau*ub;
x1=x+alpha1*s;
f1=cdmin(x1);
alpha2=tau*lb+(1-tau)*ub;
x2=x+alpha2*s;
f2=cdmin(x2);

for k=3:n
    if f1>f2
        lb=alpha1;
        alpha1=alpha2;
        f1=f2;
        alpha2=tau*lb+(1-tau)*ub;
        x2=x+alpha2*s;
        f2=cdmin(x2);
    else
        ub=alpha2;
        alpha2=alpha1;
        f2=f1;
        alpha1=(1-tau)*lb+tau*ub;
        x1=x+alpha1*s;
        f1=cdmin(x1);
    end
end    

alpha=(alpha1+alpha2)/2;
end
