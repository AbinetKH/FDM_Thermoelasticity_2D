% Title - Calculation of Heat Flow and Linear Elasticity for rectangular domain using FDM
% Written by - Abinet Kifle Habtemariam 
% Bauhaus University 
% Date 26/02/2014
clc;
clear all;
close all;
display(' _________3___________');
display('|                     |');
display('|      Boundary ID    |');
display('1(L)     for          2');
display('|     Rectangular     |');
display('|       Domain        |');
display('|                     |');
display('|__________4(W)_______|');
%  % Required input parameters for the geometrical layout 
% l= input('Enter the length for the rectangular domain [m]: ');
% w= input('Enter the width for the rectangular domain [m]: ');
% nl=input('Enter the number of mesh in l direction:');
% nw=input('Enter the number of mesh in w direction:');
% dt = input ('Time step = ');
% t_tot=input('Total time = ');
% alpha = input('thermal diffusivity [m2/s] = ');
% E = input ('Young modulus [GPa] = ');
% nu = input ('Poisson ratio= ');
% animation=input ('To create animation of the heat flow type one : ');
% scale=input('Scale factor for deformation plot =');
% display('Insert Initial condition');
% display('       a_1*y^2 + a_2*y + a_3*x^2 + a_4*x + a_5*x*y + a_6+ b_3*sin(b_1*y + b_2*x) + b_4*cos(b_5*y + b_6*x)');
%        a_1=input ('      a_1= ');a_2=input ('      a_2= ');a_3=input ('      a_3= ');
%        a_4=input ('      a_4= ');a_5=input ('      a_5= ');a_6=input ('      a_6= ');
%        b_1=input ('      b_1= ');b_2=input ('      b_2= ');b_3=input ('      b_3= ');
%        b_4=input ('      b_4= ');b_5=input ('      b_5= ');b_6=input ('      b_6= ');
% BC=zeros(1,4);c_0=zeros(1,4);c_1=zeros(1,4);c_2=zeros(1,4);c_3=zeros(1,4);c_4=zeros(1,4);
% c_5=zeros(1,4);c_6=zeros(1,4);d_1=zeros(1,4);d_2=zeros(1,4);d_3=zeros(1,4);d_4=zeros(1,4);
% d_5=zeros(1,4);d_6=zeros(1,4);d_7=zeros(1,4);d_8=zeros(1,4);d_9=zeros(1,4);d_10=zeros(1,4);d_11=zeros(1,4);
% for i=1:4;
%    display(['Insert boundary condition for side --> ' num2str(i)]);
% BC(i)=input('    Type one to choose Neumann BC or Zero for Dirichlet BC: ');
%    if BC(i)==1;
%        display('      c_0+c_1*t+c_2*t^2+c_3*sin(c_4*t)+c_5*cos(c_6*t)');
%        c_0(i)=input ('      c_0= ');c_1(i)=input ('      c_1= ');c_2(i)=input ('      c_2= ');
%        c_3(i)=input ('      c_3= ');c_4(i)=input ('      c_4= ');c_5(i)=input ('      c_5= ');
%        c_6(i)=input ('      c_6= ');
%    else
%        display('      d_1+d_2*x+d_3*x^2+d_4*cos(d_5*x+d_6*(n*dt))+d_7*sin(d_8*x+d_9*(n*dt))+d_10*(n*dt)+d_11*(n*dt)^2');
%        d_1(i)=input ('      d_1= ');d_2(i)=input ('      d_2= ');d_3(i)=input ('      d_3= ');
%        d_4(i)=input ('      d_4= ');d_5(i)=input ('      d_5= ');d_6(i)=input ('      d_6= ');
%        d_7(i)=input ('      d_7= ');d_8(i)=input ('      d_8= ');d_9(i)=input ('      d_9= ');
%        d_10(i)=input ('      d_10= ');d_11(i)=input ('      d_11= ');
%    end 
% end
%     % -------------short input-------------------------
l=3;w=5;nl=30;nw=50;alpha=4.2e-5; dt=10; t_tot=10000;E=200;nu=0.2;
% IC 
a_1=00;a_2=0;a_3=0;a_4=0;a_5=0;a_6=0;b_1=0;b_2=pi/2;b_3=0;b_4=0;b_5=pi/2;b_6=0;
% one to choose Neumann BC or Zero for Dirichlet BC
BC=[0 0 0 0];
% cofficent for Neumann BC 
c_0=[0 0 0 0];c_1=[0 0 0 0];c_2=[0 0 0 0];c_3=[0 0 0 0];c_4=[0 0 0 0];c_5=[0 0 0 0];c_6=[0 0 0 0];
% cofficent for Dirichlet BC 
d_1=[20 20 20 20];d_2=[0 0 0 0];d_3=[0 0 0 0];d_4=[0 0 0 0];d_5=[0 0 0 0];d_6=[0 0 0 0];d_7=[0 0 0 0];
d_8=[0 0 0 0];d_9=[0 0 0 0];d_10=[0 0 0 0];d_11=[0 0 0 0];
% 
animation=0;scale=500;
% BC for free(1)or fixed(2) ( for solving linear Elasticity part )
BC2=[2 2 1 1];
%      -------------- end -----------------------------------
f=zeros(1,4);
dy=l/nl;
dx=w/nw;
nt=floor(t_tot/dt);
nx=0:dy:l;
ny=0:dx:w;
size_x=size(nx,2);
size_y=size(ny,2);
figure(1)
    [X,Y] = meshgrid(ny,nx);
    surf(X,Y,zeros(size(Y)));  
    hold on
    axis equal
    view(0,90)
    xlabel('Width[m]'); ylabel('Length[m]');
    title('Meshed Rectangular Domain','FontSize',14,'Color','black','FontName', 'Times New Roman');
% solving the heat equation using implicit method
tic
T = zeros(size_x,size_y,nt);
T_new_vector = zeros(size_x*size_y,1);
A=sparse(size_x*size_y,size_x*size_y); % Sparse matrix storage , save time and storage size 
sy = alpha*dt/dy^2;
sx = alpha*dt/dx^2;
mat=zeros(size_x,size_y);
int=1;
for i=1:size_x
    for j=1:size_y
        mat(i,j) =int;
        int =int+1;
    end
end
for i=2:size_x-1
    for j=2:size_y-1
        N= mat(i,j);
        A(N,N) = 1 + 2*sx + 2*sy;
        A(N,mat(i+1,j)) = -sy;
        A(N,mat(i-1,j)) = -sy;
        A(N,mat(i,j+1)) = -sx;
        A(N,mat(i,j-1)) = -sx;   
    end
end
% boundary conditions
for j=2:size_y-1;
    A(mat(1,j),mat(1,j))=1;
    A(mat(size_x,j),mat(size_x,j))=1;
    if BC(4)==1;
        A(mat(1,j),mat(1,j))=-1;
        A(mat(1,j),mat(1,j)+size_y)=1;
    end
    if BC(3)==1;
        A(mat(size_x,j),mat(size_x,j))=1;
        A(mat(size_x,j),mat(size_x,j)-size_y)=-1;
    end  
end
for i=1:size_x;
    A(mat(i,1),mat(i,1))=1;
    A(mat(i,size_y),mat(i,size_y))=1; 
    if BC(1)==1;
        A(mat(i,1),mat(i,1))=-1;
        A(mat(i,1),mat(i,2))=1;
     end
    if BC(2)==1;
        A(mat(i,size_y),mat(i,size_y))=1;
        A(mat(i,size_y),mat(i,size_y-1))=-1;
    end
end
%Intial Condition 
y = linspace(0,l,size_x);
x = linspace(0,w,size_y);
for j=2:size_y-1
  for i=2:size_x-1
     T(i,j,1) = a_1*y(1,i)^2 + a_2*y(1,i) + a_3*x(1,j)^2 + a_4*x(1,j) + a_5*y(1,i)*x(1,j) + a_6...
     + b_2*sin(b_1*y(1,i) + b_3*x(1,j)) + b_4*cos(b_5*y(1,i) + b_6*x(1,j));
  end
end
T_new = T;
for n=2:nt
% Compute rhs
rhs = zeros(size_x*size_y,1);
for i = 2:size_x-1
   for j = 2:size_y-1
   ii = mat(i,j);
   rhs(ii)=T(i,j);
   end 
end
for i=1:4;
    if BC(i)==1
    f(1,i)=c_0(1,i)+c_1(1,i)*(n*dt)+c_2(1,i)*(n*dt)^2+c_3(1,i)*sin(c_4(1,i)*(n*dt))+c_5(1,i)*cos(c_6(1,i)*(n*dt));
    else
        f(1,i)=0;
    end
end
% at the boundaries
k_4=zeros(1,size_y-2);k_3=zeros(1,size_y-2);k_2=zeros(1,size_x-2);k_1=zeros(1,size_x-2);
for j=2:size_y-1    
 k_4(j)=d_1(4)+d_2(4)*x(1,j)+d_3(4)*x(1,j)^2+d_4(4)*cos(d_5(4)*x(1,j)+d_6(4)*(n*dt))+d_7(4)*sin(d_8(4)*x(1,j)+d_9(4)*(n*dt))+d_10(4)*(n*dt)+d_11(4)*(n*dt)^2;
 rhs(mat(1,j)) =dy*f(1,4)+k_4(j);        %boundary 4
 k_3(j)=d_1(3)+d_2(3)*x(1,j)+d_3(3)*x(1,j)^2+d_4(3)*cos(d_5(3)*x(1,j)+d_6(3)*(n*dt))+d_7(3)*sin(d_8(3)*x(1,j)+d_9(3)*(n*dt))+d_10(3)*(n*dt)+d_11(3)*(n*dt)^2;
 rhs(mat(size_x,j)) =dy*f(1,3)+k_3(j); %boundary 3
 end
for i=2:size_x-1
 k_2(i)=d_1(2)+d_2(2)*y(1,i)+d_3(2)*y(1,i)^2+d_4(2)*cos(d_5(2)*y(1,i)+d_6(2)*(n*dt))+d_7(2)*sin(d_8(2)*y(1,i)+d_9(2)*(n*dt))+d_10(2)*(n*dt)+d_11(2)*(n*dt)^2;
 rhs(mat(i,size_y))=dx*f(1,2)+k_2(i); %boundary 2
 k_1(i)=d_1(1)+d_2(1)*y(1,i)+d_3(1)*y(1,i)^2+d_4(1)*cos(d_5(1)*y(1,i)+d_6(1)*(n*dt))+d_7(1)*sin(d_8(1)*y(1,i)+d_9(1)*(n*dt))+d_10(1)*(n*dt)+d_11(1)*(n*dt)^2;
 rhs(mat(i,1)) =dx*f(1,1)+k_1(i);  %boundary 1
end  
% at the boundaries
% rhs(mat(1,1:size_y)) =dy*f(1,4)+DBC(1,4);        %boundary 4
% rhs(mat(size_x,1:size_y)) =dy*f(1,3)+DBC(1,3);   %boundary 3
% rhs(mat(2:size_x-1,1)) =dx*f(1,1)+DBC(1,1);       %boundary 1
% rhs(mat(2:size_x-1,size_y))=dx*f(1,2)+DBC(1,2);   %boundary 2

T_new_vector=A\rhs;
T_new(:,:,n) = T_new_vector(mat);
T = T_new(:,:,n);
% corner values 
T(1,1) = (T(2,1)+T(1,2))/2;
T(size_x,1) = (T(size_x-1,1)+T(size_x,2))/2;
T(1,size_y) = (T(2,size_y)+T(1,size_y-1))/2;
T(size_x,size_y) = (T(size_x-1,size_y)+ T(size_x,size_y-1))/2;

if animation==1
   if (mod(n,20)==0)
       set(gcf,'color', 'white', 'Units', 'normalized','Position', [0.25, 0.25, 0.50, 0.50])
     clf    
    surf(X,Y,T);
    %axis([0 w 0 l 0 max(T_new_vector)]);
    %view(0,90);
    xlabel('width [m]'); ylabel('length[m]'); zlabel('Temperature [°]');
    drawnow update
    format short
    title({'Heat Flow-Implicit Method';[' Time ',num2str(n*dt),'s','(',num2str(round((n*dt)/60)),'min)']},'FontSize',16,'Color','blue','FontName', 'Times New Roman');
        %legend(['        Time ',num2str(n*dt),'s'],'Location','NorthWest')
    colorbar('location','eastoutside','fontsize',12);
    F=getframe;
   end
end
end
% movie(F,fram,1)
figure(2)
    set(gcf, 'Units', 'normalized', 'Position', [0.25, 0.25, 0.50, 0.50])
    clf
    surfc(X,Y,T);
    %shading interp;
    %imagesc(T);
%     axis('equal')
    xlabel('Length [m]'); ylabel('Height [m]'); zlabel('Temperature [°]');
    title('Heat Flow-Implicit Method','FontSize',16,'Color','blue','FontName', 'Times New Roman');
    drawnow
    colorbar

toc
%% linear Elasticity Equation
% solve disp. B*u=rhs
%BC2=[2 2 1 1]; %short input 
% BC2=zeros(1,4);
% display('Specify the boundary conditon for solving linear Elasticity Equation');
%   display('  1 for free BC'); 
%   display('  2 for fixed BC'); 
% for i=1:4
%   display(['Insert boundary condition for side --> ' num2str(i)]);
%      BC2(i)=input('    Choose one of the boundary condition by typing the number: ');  
%     if BC2(i)==1;
%         display(['      You have chosen free boundary condition for side--> ' num2str(i)]);
%     end
%     if BC2(i)==2;
%         display(['      You have chosen fixed boundary condition for side--> ' num2str(i)]);
%     end
% end
display('The computation might take some minutes... ');
tic
% Lamé parameters
lamda=(E*nu)/((1+nu)*(1-2*nu));
miu=E/(2*(1+nu));

% Construct the B matrix 
B = sparse(size_x*size_y*2,size_x*size_y*2);%sparse
b = lamda + 2*miu;
d = lamda + miu;
c = miu;
s=size_y*size_x;
for i=2:size_x-1
    for j=2:size_y-1
        N2= mat(i,j);
        %x direction(U) 
        B(N2,N2)=-((2*b/dx^2)+(2*c/dy^2));
        B(N2,mat(i-1,j))=c/dy^2;
        B(N2,mat(i+1,j))=c/dy^2;
        B(N2,mat(i,j-1))=b/dx^2 ;
        B(N2,mat(i,j+1))=b/dx^2 ;
        B(N2,mat(i+1,j+1)+s)=d/(4*dx*dy);
        B(N2,mat(i+1,j-1)+s)=-d/(4*dx*dy);
        B(N2,mat(i-1,j+1)+s)=-d/(4*dx*dy);
        B(N2,mat(i-1,j-1)+s)=d/(4*dx*dy);
               
        %y direction(V)
        B(N2+s,N2+s)=-((2*b/dy^2)+(2*c/dx^2));
        B(N2+s,mat(i,j-1)+s)=c/dx^2;
        B(N2+s,mat(i,j+1)+s)=c/dx^2;
        B(N2+s,mat(i-1,j)+s)=b/dy^2;
        B(N2+s,mat(i+1,j)+s)=b/dy^2;
        B(N2+s,mat(i+1,j+1))=d/(4*dx*dy);
        B(N2+s,mat(i+1,j-1))=-d/(4*dx*dy);
        B(N2+s,mat(i-1,j+1))=-d/(4*dx*dy);
        B(N2+s,mat(i-1,j-1))=d/(4*dx*dy);
      
    end
end

% Boundary condition 
   
d1=lamda;
d3=E/(2*(1+nu));

for i=2:size_x-1;
  if BC2(1,1)==1;
    % side 1 in u or x direction 
    B(mat(i,1),mat(i+1,1)+s)=d1/(2*dy);
    B(mat(i,1),mat(i-1,1)+s)=-d1/(2*dy);
    B(mat(i,1),mat(i,1))=-b/dx;
    B(mat(i,1),mat(i,2))=b/dx;
    % side 1 in v or y direction 
    B(mat(i,1)+s,mat(i,2)+s)=d3/dx;
    B(mat(i,1)+s,mat(i,1)+s)=-d3/dx;
    B(mat(i,1)+s,mat(i-1,1))=-d3/(2*dy);
    B(mat(i,1)+s,mat(i+1,1))=d3/(2*dy);
  end
  if BC2(1,2)==1;
    % side 2 in u or x direction 
    B(mat(i,size_y),mat(i-1,size_y)+s)=-d1/(2*dy);
    B(mat(i,size_y),mat(i+1,size_y)+s)=d1/(2*dy);
    B(mat(i,size_y),mat(i,size_y))=b/dx;
    B(mat(i,size_y),mat(i,size_y-1))=-b/dx;
    % side 2 in v or y direction 
    B(mat(i,size_y)+s,mat(i,size_y-1)+s)=-d3/dx;
    B(mat(i,size_y)+s,mat(i,size_y)+s)=d3/dx;
    B(mat(i,size_y)+s,mat(i+1,size_y))=d3/(2*dy);
    B(mat(i,size_y)+s,mat(i-1,size_y))=-d3/(2*dy);
   end      
end
for j=2:size_y-1;
    if BC2(1,4)==1;
    % side 4 in v or y direction 
    B(mat(1,j)+s,mat(1,j+1))=d1/(2*dx);
    B(mat(1,j)+s,mat(1,j-1))=-d1/(2*dx);
    B(mat(1,j)+s,mat(1,j)+s)=-b/dy;
    B(mat(1,j)+s,mat(2,j)+s)=b/dy;
    % side 4 in u or x direction 
    B(mat(1,j),mat(1,j+1)+s)=d3/(2*dx);
    B(mat(1,j),mat(1,j-1)+s)=-d3/(2*dx);
    B(mat(1,j),mat(2,j))=d3/dy;
    B(mat(1,j),mat(1,j))=-d3/dy;
    end 
    if BC2(1,3)==1;
    % side 3 in v or y direction 
    B(mat(size_x,j)+s,mat(size_x,j+1))=d1/(2*dx);
    B(mat(size_x,j)+s,mat(size_x,j-1))=-d1/(2*dx);
    B(mat(size_x,j)+s,mat(size_x,j)+s)=b/dy;
    B(mat(size_x,j)+s,mat(size_x-1,j)+s)=-b/dy;
    % side 3 in u or x direction 
    B(mat(size_x,j),mat(size_x,j+1)+s)=d3/(2*dx);
    B(mat(size_x,j),mat(size_x,j-1)+s)=-d3/(2*dx);
    B(mat(size_x,j),mat(size_x-1,j))=-d3/dy;
    B(mat(size_x,j),mat(size_x,j))=d3/dy;
    end

end

%BC for corner points
if BC2(1,4)==1==BC2(1,2);
% side 4 right node in u direction 
B(mat(1,size_y),mat(1,size_y))=b/dx+d3/dy;
B(mat(1,size_y),mat(1,size_y)+s)=-d1/dy-d3/dx;
B(mat(1,size_y),mat(1,size_y-1))=-b/dx;
B(mat(1,size_y),mat(2,size_y)+s)=d1/dy;
B(mat(1,size_y),mat(2,size_y))=-d3/dy;
B(mat(1,size_y),mat(1,size_y-1)+s)=d3/dx;
end
if BC2(1,4)==1==BC2(1,2);
% side 4 right node in v direction 
B(mat(1,size_y)+s,mat(1,size_y))=-d1/dx-d3/dy;
B(mat(1,size_y)+s,mat(1,size_y)+s)=d3/dx+b/dy;
B(mat(1,size_y)+s,mat(1,size_y-1))=d1/dx;
B(mat(1,size_y)+s,mat(2,size_y)+s)=-b/dy;
B(mat(1,size_y)+s,mat(2,size_y))=d3/dy;
B(mat(1,size_y)+s,mat(1,size_y-1)+s)=-d3/dx;
end
if BC2(1,3)==1==BC2(1,2);
% side 3 right node in u direction 
B(mat(size_x,size_y),mat(size_x,size_y))=b/dx+d3/dy;
B(mat(size_x,size_y),mat(size_x,size_y)+s)=d1/dy+d3/dx;
B(mat(size_x,size_y),mat(size_x,size_y-1))=-b/dx;
B(mat(size_x,size_y),mat(size_x-1,size_y)+s)=-d1/dy;
B(mat(size_x,size_y),mat(size_x-1,size_y))=-d3/dy;
B(mat(size_x,size_y),mat(size_x,size_y-1)+s)=-d3/dx;
end
if BC2(1,3)==1==BC2(1,2);
% side 3 right node in v direction 
B(mat(size_x,size_y)+s,mat(size_x,size_y))=d1/dx+d3/dy;
B(mat(size_x,size_y)+s,mat(size_x,size_y)+s)=d3/dx+b/dy;
B(mat(size_x,size_y)+s,mat(size_x,size_y-1))=-d1/dx;
B(mat(size_x,size_y)+s,mat(size_x-1,size_y)+s)=-b/dy;
B(mat(size_x,size_y)+s,mat(size_x-1,size_y))=-d3/dy;
B(mat(size_x,size_y)+s,mat(size_x,size_y-1)+s)=-d3/dx;
end
if BC2(1,4)==1==BC2(1,1);
% side 4 left node in u direction -
B(mat(1,1),mat(1,1))=b/dx+d3/dy;
B(mat(1,1),mat(1,1)+s)=d1/dy+d3/dx;
B(mat(1,1),mat(1,2))=-b/dx;
B(mat(1,1),mat(2,1)+s)=-d1/dy;
B(mat(1,1),mat(2,1))=-d3/dy;
B(mat(1,1),mat(1,2)+s)=-d3/dx;
end
if BC2(1,4)==1==BC2(1,1);
% side 4 left node in v direction -
B(mat(1,1)+s,mat(1,1))=d1/dx+d3/dy;
B(mat(1,1)+s,mat(1,1)+s)=d3/dx+b/dy;
B(mat(1,1)+s,mat(1,2))=-d1/dx;
B(mat(1,1)+s,mat(2,1)+s)=-b/dy;
B(mat(1,1)+s,mat(2,1))=-d3/dy;
B(mat(1,1)+s,mat(1,2)+s)=-d3/dx;
end
if BC2(1,3)==1==BC2(1,1);
% side 3 left node in u direction 
B(mat(size_x,1),mat(size_x,1))=b/dx+d3/dy;
B(mat(size_x,1),mat(size_x,1)+s)=-d1/dy-d3/dx;
B(mat(size_x,1),mat(size_x,2))=-b/dx;
B(mat(size_x,1),mat(size_x-1,1)+s)=d1/dy;
B(mat(size_x,1),mat(size_x-1,1))=-d3/dy;
B(mat(size_x,1),mat(size_x,2)+s)=d3/dx;
end
if BC2(1,3)==1==BC2(1,1);
% side 3 left node in v direction 
B(mat(size_x,1)+s,mat(size_x,1))=-d1/dx-d3/dy;
B(mat(size_x,1)+s,mat(size_x,1)+s)=d3/dx+b/dy;
B(mat(size_x,1)+s,mat(size_x,2))=d1/dx;
B(mat(size_x,1)+s,mat(size_x-1,1)+s)=-b/dy;
B(mat(size_x,1)+s,mat(size_x-1,1))=d3/dy;
B(mat(size_x,1)+s,mat(size_x,2)+s)=-d3/dx;

end

% Displacment BC 
if BC2(1,1)==2
   for k=1:size_x
    B(mat(k,1),mat(k,1))=1;
    B(mat(k,1)+s,mat(k,1)+s)=1;
   end
end 
if BC2(1,2)==2
   for k=1:size_x
    B(mat(k,size_y),mat(k,size_y))=1;
    B(mat(k,size_y)+s,mat(k,size_y)+s)=1;
   end
end 
if BC2(1,3)==2
   for k=1:size_y
    B(mat(size_x,k),mat(size_x,k))=1;
    B(mat(size_x,k)+s,mat(size_x,k)+s)=1;
   end
end 
if BC2(1,4)==2
   for k=1:size_y
    B(mat(1,k),mat(1,k))=1;
    B(mat(1,k)+s,mat(1,k)+s)=1;
   end
end 
alfa2=1;
% Compute rhs2
rhs2 = zeros((2*size_x*size_y),1);
    for i = 2:size_x-1
        for j = 2:size_y-1
%             rhs2(mat(i,j)) = *((T_new(i,j+1,nt)-T_new(i,j-1,nt))/(2*dy)-(T_new(i,j+1,1)-T_new(i,j-1,1))/(2*dy));
%             rhs2(mat(i,j)+s) = ((T_new(i+1,j,nt)-T_new(i-1,j,nt))/(2*dx)-(T_new(i,j+1,1)-T_new(i,j-1,1))/(2*dx));
           % rhs2(mat(i,j)) = alpha*(3*lamda+2*miu)*((T_new(i,j+1,nt)-T_new(i,j-1,nt))/(2*dx));%-(T_new(i,j+1,1)-T_new(i,j-1,1))/(2*dx));
            %rhs2(mat(i,j)+s) = alpha*(3*lamda+2*miu)*((T_new(i+1,j,nt)-T_new(i-1,j,nt))/(2*dy));%-(T_new(i,j+1,1)-T_new(i,j-1,1))/(2*dy));
            rhs2(mat(i,j)) = alpha*(3*lamda+2*miu)*((T_new(i,j+1,nt)-T_new(i,j-1,nt))/(2*dx)-(T_new(i,j+1,1)-T_new(i,j-1,1))/(2*dx));
            rhs2(mat(i,j)+s) = alpha*(3*lamda+2*miu)*((T_new(i+1,j,nt)-T_new(i-1,j,nt))/(2*dy)-(T_new(i,j+1,1)-T_new(i,j-1,1))/(2*dy));
        end
    end
%   for i=2:size_x-1
%       rhs2(mat(i,size_y))=((T_new(i,size_y,nt)-T_new(i,size_y-1,nt))/(2*dx)-(T_new(i,size_y,1)-T_new(i,size_y-1,1))/(2*dx));
%       rhs2(mat(i,1))=((-T_new(i,1,nt)+T_new(i,2,nt))/(2*dx)-(-T_new(i,1,1)+T_new(i,2,1))/(2*dx));
%   end
%   for j=1:size_y
%      rhs2(mat(size_x,j))=((T_new(size_x,j,nt)-T_new(size_x-1,j,nt))/(2*dy)-(T_new(size_x,j,1)-T_new(size_x-1,j,1))/(2*dy)); 
%      rhs2(mat(1,j))=((-T_new(1,j,nt)+T_new(2,j,nt))/(2*dy)-(-T_new(1,j,1)+T_new(2,j,1))/(2*dy)); 
%   end
% Compute solution vector
UV= B\rhs2;

% Create the matrix with displacements for both directions
U_m = UV(mat(1:size_x,1:size_y));
V_m = UV(mat(1:size_x,1:size_y)+s);

figure(3)
hold on
quiver(U_m,V_m);
title('Gradient of Displacement','FontSize',16,'Color','black','FontName','Times New Roman');
xlabel('Number of nodes x'); ylabel('Number of nodes y');
%Deformed shape
U = UV(1:s,1);
V = UV(s+1:2*s,1);
U_mm=vec2mat(U,size_y);
V_mm=vec2mat(V,size_y);
for i=1:size_x
    for j=1:size_y
        Xd(i,j)=X(i,j)+U_mm(i,j)*scale;
        Yd(i,j)=Y(i,j)+V_mm(i,j)*scale;
     end
end 
figure(4)
   surf(Xd,Yd,zeros(size(Y)));  
   view(0,90);
   hold on
   rectangle('Position',[1e-5,1e-5,w,l],'EdgeColor','r','LineWidth',2,'LineStyle','--')
   xlabel('Width[m]'); ylabel('Length[m]');
   title('Deformed Shape','FontSize',16,'Color','black','FontName', 'Times New Roman');
toc








































