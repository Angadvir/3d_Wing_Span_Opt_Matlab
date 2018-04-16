%% Forward Difference Gradient Approximator


function [gradf]=fwd(fcn,x,r)

pf=fcn(x);

for i=1:length(x)

    dx=abs(r*x(i));
    
    xtemp1=x;
    xtemp1(i)=xtemp1(i)+dx;

    pfpx(i)=(fcn(xtemp1)-pf)/(dx);
    
end

gradf=pfpx;

end


