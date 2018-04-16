function[gSup, hSup] = Wingdrag_Con(x)
global tau, global tcap, global rh, global g, global t, global Mr, global wWeb
global Nlift,  global Lb, global wbar,global sigmaMaxSh,global Mrb, global wCap
global tweb,  global sigmaMax, global rhoWeb, global rhoCap, global Icap

gSup = [((Lb*t*x(1)*(x(4)^2/x(5))/(24*x(2)))-Mr);                                       %g1
        (Nlift*(x(4)^2)* Mrb*tau*x(2)^2/(x(5)^2)*(Icap*sigmaMax))-8;                    %g3 
        (((x(4)^2*Lb*Nlift*x(2)^2)/(tau*(x(5)^2)*tweb*sigmaMaxSh)-12));                 %g4 
        (0.86*(x(1)^-2.38) + 0.14*(x(1)^0.56)-(x(6)^3.94));                   %Auxiliarry Constraint
        ((8*rhoCap*g*wbar*tcap*(x(5)^1.5)*x(6))/(3*(x(4)^2/x(5))^0.5)-wCap);            %g5
        (x(6)*(8*rhoWeb*g*rh*tau*tweb*(x(5)^1.5)*x(6))/(3*(x(4)^2/x(5))^0.5)-wWeb)];    %g6

hSup =  [(-x(5) + (x(4)/4)*((t/(x(3)*tau))+(t/tau)));                      %h1
         ((x(4)^2/x(5))-(x(4)^2/((x(4)/4)*((t/(x(3)*tau))+(t/tau)))));     %h2    
          x(1)-1-(2*x(3));                                                 %h3    
          x(2)-1-x(3);                                                     %h4
          x(6)-(x(3)^2+x(3)+1)/(1+x(3))^2];                                %h5
end 
%%
%X1[x(1),x(2),x(3),x(4),x(5),x(6)]=[p,q,Lambda,b,v(Lambda)]