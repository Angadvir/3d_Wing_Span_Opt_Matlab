%% Conjugate Gradient Optimizer

function [x]=conjugate(fcn,x0,grad_method,r, n,epsilon, nmax)
tic

% For the optimization of a NACA 0012 airfoil. 
fcn=@cdmin;
x0=[0.0146,0.3025,0.06,-0.4928,0.3016,0.06,-0.4848,-0.0039,0.0101,-2.7791,9.2496];

grad_method=1;
r=1e-3;
n=27;
epsilon=1e-8;
nmax=100;

% Initial x
x=x0;
i=0;
ii=0;

disp('Fletcher-Reeves Conjugate Gradient Method')

F(1)=fcn(x);
formatSpec='i: %d    ii=%d    f(x): %.6e \n';
fprintf(formatSpec,i,ii,F(1));%,xtr);
figure('color','white')

% Initial Gradient 
        if grad_method==1
            gradf=central(fcn,x,r);
        elseif grad_method==2
            gradf=fwd(fcn,x,r);
        else
            error('Error: select gradient approximation method.')
        end

% initialize
a=dot(gradf,gradf);     % Dot product

% Loop
i=0;    % iteration counter
ii=0;    % Counter for inner while loop

%alpha=1;
for i=1:nmax
    s=-gradf;
    alpha=golden(x,gradf,n);
    
    if abs(alpha)<1e-6  
        disp('Optimization terminated: alpha=0')
        break
    else
        x=x+alpha.*s;
    end
    
    while max(abs(gradf))>epsilon
        if grad_method==1
            gradf=central(fcn,x,r);
        elseif grad_method==2
            gradf=fwd(fcn,x,r);
        else
            error('Error: select gradient approximation method.')
            break
        end
    b=dot(gradf,gradf);
    beta=b/a;
    s=-gradf+beta*s;
    a=b;
    slope=dot(s,gradf);

    if slope>=0
        %disp('Break: Slope')
        break                   % Breaks while loop
    else
        alpha=golden(x,gradf,n);
        x=x+alpha.*s;
        
        %% Display: Comment out for improved performance
        ii=ii+1;
        F(ii+1)=fcn(x);
        formatSpec='i: %d    ii=%d    f(x): %.6e \n';
        fprintf(formatSpec,i,ii,F(ii+1));
        
    end
    end
            ii=ii+1;
            F(ii+1)=fcn(x);
    if max(abs(gradf))<epsilon   %Break for loop
        
        f=fcn(x);      
        formatSpec='Converged:abs(gradient)<epsilon\nf(x): %.6e\n\n';
        fprintf(formatSpec,f);
        if F==fcn(x0)
        disp('***Initial point is a local minimum.***')
        disp('Check n in golden.m and dx in central.')
        end
        break
    end
    if i==nmax
        disp('Number of function evaluations exceeded.')
    end
end
toc
end

