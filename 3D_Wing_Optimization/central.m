%% Central Difference Gradient Approximation


function [gradf]=central(fcn,x,r)



for i=1:length(x)
    
    dx=abs(r*x(i));
    
    xtemp1=x;
    xtemp2=x;
    xtemp1(i)=xtemp1(i)+dx;
    xtemp2(i)=xtemp2(i)-dx;

    pfpx(i)=(fcn(xtemp1)-fcn(xtemp2))/(2*dx);
    
end

gradf=pfpx;

end

