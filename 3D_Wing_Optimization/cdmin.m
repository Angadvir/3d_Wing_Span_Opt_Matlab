%% Airfoil Optimization Objective Function
function cd = cdmin(p)

a=parsec(p);

point=100;
n=(1/point);
nn=n/5; % to create more points at the LE

% upper surface coordinates
x1=1:-n:0.1;
x2=0.1:-nn:0;

% lower surface coordinates
%xxx=0:n:1;
x3=0:nn:0.1;
x4=0.1:n:1;

xu(1:length(x1))=x1;
xu(length(x1)+1:length(x1)+length(x2))=x2;

xl(1:length(x3))=x3;
xl(length(x3)+1:length(x3)+length(x4))=x4;

X(1:length(xu))=xu;
X(length(xu)+1:length(xu)+length(xl))=xl;

X=X';

yu=a(1)*xu.^.5+a(2)*xu.^(1.5)+a(3)*xu.^(2.5)+a(4)*xu.^(3.5)+a(5)*xu.^(4.5)+a(6)*xu.^(5.5);
yl=-(a(7)*xl.^.5+a(8)*xl.^(1.5)+a(9)*xl.^(2.5)+a(10)*xl.^(3.5)+a(11)*xl.^(4.5)+a(12)*xl.^(5.5));

Y(1:length(yu),1)=yu;
Y(length(yu)+1:length(yu)+length(yl),1)=yl;

Y=round(Y,7);
X=round(X,7);
coord=[X Y];

% Run Xfoil, requires xfoil.m
[a,b]=xfoil(coord,5E5,0.5);

% Plot airfoil evolution, comment out if desired
%hold on
plot(X,Y)
%figure
grid minor
hold on
grid minor
axis equal
        
% Objective function
cd=a.CD; 
end



