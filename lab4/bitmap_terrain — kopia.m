% This is a very fast implementation of bitmap terrain projection with direction splitting.
% It caches basis functions values for integration points prior to main integration loop.
%
% How to use
% 
% bitmap_terrain(filename as a string, number of elements along x axis, polynomial order alog x axis, number of elements along y axis, polynomial order alog y axis)
%
% Examples
% 
% bitmap_terrain(129,"C:\\Users\\mpasz\\Documents\\Terrain.png",62,2,62,2)
% The precision here is 2(nx+2)+1=2*(62+2)+1 which gives two points per element in one direction for a final plot

%In this updated version we only provide file location, 
% number of elements in each direction, and the B-spline order
function [GRAY] = bitmap_terrain_error(filename,elements,p)
  
tic;
elementsx=elements;
elementsy=elements;
px=p; py=p;
precision = 2*(elements+p)+1;
% subroutine calculating number of basis functions
compute_nr_basis_functions = @(knot_vector,p) size(knot_vector, 2) - p - 1;
% subroutine calculating mesh for plotting splines
mesh = @(a,c,precision) [a:(c-a)/precision:c];

% create knot_vectors along x any y axis
knot_vectorx = simple_knot(elementsx,px);
knot_vectory = simple_knot(elementsy,py);

% read image
X = imread(filename); %to make it upside down just add 255 - 
X = flipud(X);
% extract R, G and B components of the image
R = 255-X(:,:,1);
G = 255-X(:,:,2);
B = 255-X(:,:,3);

% read size of image
ix = size(X,1);
iy = size(X,2);

% compute number of degrees of freedom
nx = number_of_dofs(knot_vectorx,px);
ny = number_of_dofs(knot_vectory,py);

% initiate matrices for further computations
Ax = sparse(nx,nx);
Ay = sparse(ny,ny);
FRx = zeros(nx,ny);
FGx = zeros(nx,ny);
FBx = zeros(nx,ny);

init = toc
tic;

% initiate matrices for precached basis function values at given points
splinex = zeros(elementsx,nx,px+1);
spliney = zeros(elementsy,ny,py+1);

% precache values of basis functions in integration points

for ex = 1:elementsx;
% range of nonzero functions over element
  [xl,xh] = dofs_on_element(knot_vectorx,px,ex);
% range of element (left and right edge over x axis)
  [ex_bound_l,ex_bound_h] = element_boundary(knot_vectorx,px,ex);
% quadrature points over element (over x axis)
  qpx = quad_points(ex_bound_l,ex_bound_h,px+1);
% quadrature weights over element (over x axis)
  qwx = quad_weights(ex_bound_l,ex_bound_h,px+1);
% loop over nonzero functions over element
  for bi = xl:xh
% loop over quadrature points
    for iqx = 1:size(qpx,2)
      splinex(ex,bi,iqx)= compute_spline(knot_vectorx,px,bi,qpx(iqx));
    end
  end
end

for ey = 1:elementsy;
% range of nonzero functions over element
  [yl,yh] = dofs_on_element(knot_vectory,py,ey);
% range of element (left and right edge over y axis)
  [ey_bound_l,ey_bound_h] = element_boundary(knot_vectory,py,ey);
% quadrature points over element (over y axis)
  qpy = quad_points(ey_bound_l,ey_bound_h,py+1);
% quadrature weights over element (over y axis)
  qwy = quad_weights(ey_bound_l,ey_bound_h,py+1);
% loop over nonzero functions over element
  for bi = yl:yh
% loop over quadrature points
    for iqy = 1:size(qpy,2)
      spliney(ey,bi,iqy)= compute_spline(knot_vectory,py,bi,qpy(iqy));
    end
  end
end

init_splines=toc
tic;

% integral B^x_i(x) B^y_j(y) B^x_k(x) B^y_l(y)
% (i,k=1,...,Nx; j,l=1,...,Ny)px
% loop over elements over x axis
for ex = 1:elementsx;
% range of nonzero functions over element
  [xl,xh] = dofs_on_element(knot_vectorx,px,ex);
% range of element (left and right edge over x axis)
  [ex_bound_l,ex_bound_h] = element_boundary(knot_vectorx,px,ex);
% Jacobian = size of element
  J = ex_bound_h - ex_bound_l;
% quadrature points over element (over x axis)
  qpx = quad_points(ex_bound_l,ex_bound_h,px+1);
% quadrature weights over element (over x axis)
  qwx = quad_weights(ex_bound_l,ex_bound_h,px+1);
% loop over nonzero functions over element
  for bi = xl:xh
    for bk = xl:xh
% loop over quadrature points
      for iqx = 1:size(qpx,2)
% B^x_k(x)
        funk = splinex(ex,bk,iqx);
% B^x_i(x)
        funi = splinex(ex,bi,iqx);
% B^x_i(x) B^y_j(y) B^x_k(x) B^y_l(y)
        fun = funi*funk;
        
% integral z B^x_i(x) B^y_j(y) B^x_k(x) B^y_l(y)
% (i,k=1,...,Nx; j,l=1,...,Ny)
        int = fun*qwx(iqx)*J;
        if (int~=0)
          Ax(bi,bk) = Ax(bi,bk) + int;
        end
      end
    end
  end
end

lhsx=toc


% integral B^x_i(x) B^y_j(y) B^x_k(x) B^y_l(y)
% (i,k=1,...,Nx; j,l=1,...,Ny)
% loop over elements on y axis
for ey = 1:elementsy
% range of nonzero functions over element
  [yl,yh] = dofs_on_element(knot_vectory,py,ey);
% range of element (left and right edge over y axis)
  [ey_bound_l,ey_bound_h] = element_boundary(knot_vectory,py,ey);
% Jacobian = size of element
  J = ey_bound_h - ey_bound_l;
% quadrature points over element (over y axis)
  qpy = quad_points(ey_bound_l,ey_bound_h,py+1);
% quadrature weights over element (over y axis)
  qwy = quad_weights(ey_bound_l,ey_bound_h,py+1);
% loop over nonzero functions over element
  for bj = yl:yh
    for bl = yl:yh
% loop over quadrature points      
      for iqy = 1:size(qpy,2)
% B^y_l(y)
        funl = spliney(ey,bl,iqy);
% B^y_j(y)
        funj = spliney(ey,bj,iqy);
% B^x_i(x) B^y_j(y) B^x_k(x) B^y_l(y)
        fun = funj*funl;
        
% integral z B^x_i(x) B^y_j(y) B^x_k(x) B^y_l(y)
% (i,k=1,...,Nx; j,l=1,...,Ny)
        int = fun*qwy(iqy)*J;
        if (int~=0)
          Ay(bj,bl) = Ay(bj,bl) + int;
        end
      end
    end
  end
end


lhsy=toc
tic;

            
% Integral BITMAP(x,y) B^x_k(x) B^y_l(y)
% loop over elements on x axis
for ex = 1:elementsx;
% range of nonzero functions over element
  [xl,xh] = dofs_on_element(knot_vectorx,px,ex);
% range of element (left and right edge over x axis)
  [ex_bound_l,ex_bound_h] = element_boundary(knot_vectorx,px,ex);
% loop over elements on y axis
  for ey = 1:elementsy
% range of nonzero functions over element
    [yl,yh] = dofs_on_element(knot_vectory,py,ey);
% range of element (left and right edge over y axis)
    [ey_bound_l,ey_bound_h] = element_boundary(knot_vectory,py,ey);
% Jacobian = size of element
    Jx = ex_bound_h - ex_bound_l;
    Jy = ey_bound_h - ey_bound_l;
    J = Jx * Jy;
% quadrature points over element (over x axis)
    qpx = quad_points(ex_bound_l,ex_bound_h,px+1);
% quadrature points over element (over y axis)
    qpy = quad_points(ey_bound_l,ey_bound_h,py+1);
% quadrature weights over element (over x axis)
    qwx = quad_weights(ex_bound_l,ex_bound_h,px+1);
% quadrature weights over element (over y axis)
    qwy = quad_weights(ey_bound_l,ey_bound_h,py+1);
% loop over nonzero functions over element
    for bk = xl:xh
      for bl = yl:yh  
% loop over quadrature points      
        for iqx = 1:size(qpx,2)
          for iqy = 1:size(qpy,2)
% B^x_k(x)
            funk = splinex(ex,bk,iqx);
% B^y_l(y)
            funl = spliney(ey,bl,iqy);
% integral BITMAP(x,y) B^x_k(x) B^y_l(y) over RGB components
            intR = funk*funl*qwx(iqx)*qwy(iqy)*J*bitmp(R,qpx(iqx),qpy(iqy));
            intG = funk*funl*qwx(iqx)*qwy(iqy)*J*bitmp(G,qpx(iqx),qpy(iqy));
            intB = funk*funl*qwx(iqx)*qwy(iqy)*J*bitmp(B,qpx(iqx),qpy(iqy));
            FRx(bk,bl) = FRx(bk,bl) + intR;
            FGx(bk,bl) = FGx(bk,bl) + intG;
            FBx(bk,bl) = FBx(bk,bl) + intB;
          end
        end
      end
    end
  end
end

rhs=toc
tic;

% solve one direction
[RRx,GGx,BBx]=solve_direction(Ax,FRx,FGx,FBx);

factorx=toc
tic;

% transpose matrices to solve over the other direction
FRy = transpose(RRx);
FGy = transpose(GGx);
FBy = transpose(BBx);

reorder=toc
tic

% solve second direction
[RRy,GGy,BBy]=solve_direction(Ay,FRy,FGy,FBy);

factory=toc
tic;

% transpose matrices back
RR = transpose(RRy);
GG = transpose(GGy);
BB = transpose(BBy);

% reconstruction of image

% set zero to reconstructed image matrices
R1 = zeros(ix,iy);
G1 = zeros(ix,iy);
B1 = zeros(ix,iy);


funx_tab = zeros(nx,ix);
funy_tab = zeros(ny,iy);

% precache basis functions values

% loop over basis functions
for bi = 1:nx
% loop over nonzero pixels over given function
  for i=xx(knot_vectorx(bi)):xx(knot_vectorx(bi+px+1))
% scale coordinates [1-width] -> [0-1]
    ii = (i-1)/(ix-1);
% B^x_i(x)
    funx_tab(bi,i) = compute_spline(knot_vectorx,px,bi,ii);
  end
end

% loop over basis functions
  for bj = 1:ny
% loop over nonzero pixels over given function
    for j=yy(knot_vectory(bj)):yy(knot_vectory(bj+py+1))  
% scale coordinates [1-height] -> [0-1]
      jj = (j-1)/(iy-1);
% B^y_j(y)
      funy_tab(bj,j) = compute_spline(knot_vectory,py,bj,jj);
  end
end

preprocess=toc
tic;

% terrain reconstruction
x_begin = knot_vectorx(1);
y_begin = knot_vectory(1);
% end of drawing range
x_end = knot_vectorx(size(knot_vectorx,2));
y_end = knot_vectory(size(knot_vectory,2));
x=mesh(x_begin,x_end,precision);
y=mesh(y_begin,y_end,precision);

%X and Y coordinates of points over the 2D mesh
[X,Y]=meshgrid(x,y);

%RR is 1:nx, GRAY is 1:precision
GRAY=zeros(precision+1,precision+1);
if nx==precision
  for i=1:precision
    for j=1:precision
      GRAY(i,j)=255.0-(0.3*RR(i,j)+0.59*GG(i,j)+0.11*BB(i,j));
    end
  end
else
  for i=1:precision
    for j=1:precision
      i1 = floor((i-1)/(floor(precision/nx)))+1;
      if(i1>nx)
        i1=nx;
      end
      j1 = floor((j-1)/(floor(precision/ny)))+1;
      if(j1>ny)
        j1=ny;
      end
      GRAY(i,j)=255.0-(0.3*RR(i1,j1)+0.59*GG(i1,j1)+0.11*BB(i1,j1));
    end
  end
end
for i=1:precision+1
  GRAY(precision+1,i)=GRAY(precision,i);
  GRAY(i,precision+1)=GRAY(i,precision);
end

%[U,S,V]=svd(GRAY);
%count=0;
%for i=1:precision+1
%  abs(S(i,i))
%  if(abs(S(i,i))<10)
%    S(i,i)=0.0;
%    count=count+1;
%  end
%end
%disp('Liczba sigm');
%precision+1
%disp('Niezerowe sigmy');
%precision-count+1

%GRAY = U*S*V';


% MY CHANGES
hold on;
Z = zeros(precision+1, precision+1);
for i = 1:nx
    vx = compute_spline(knot_vectorx, px, i, X);
    for j = 1:ny
        vy = compute_spline(knot_vectory, py, j, Y);
        Z = Z + vx .* vy .* GRAY;
    end
end

% Rysowanie wykresu
surf(X, Y, Z);

% Zapis rzutu z góry
view(2);
saveas(gcf, 'top_view.png');

% Zapis kilku rzutów perspektywicznych
angles = [30, 45, 60];
for k = 1:length(angles)
    view(angles(k), 30);
    saveas(gcf, sprintf('perspective_view_%d.png', angles(k)));
end


% Tworzenie animacji 360-stopniowego obrotu
frames = 360;
im = cell(1, frames);
for idx = 1:frames
    theta = (idx - 1) * (360 / frames);
    view(theta, 30);
    drawnow;
    frame = getframe(gcf);
    im{idx} = frame2im(frame);
end

% Zapis animacji do GIFa
filename = 'rotation_animation.gif';
for idx = 1:frames
    [A, map] = rgb2ind(im{idx}, 256);
    if idx == 1
        imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.05);
    else
        imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
    end
end


hold off;

% END MY CHANGES


rebuild=toc;


% Subroutine to solve one direction as 1D problem with multiple RHS
function [RR,GG,BB]=solve_direction(A,FR,FG,FB)
% compute LU factorization of A  matrix
  [L,U,P,Q] = lu(A);
  Q1=Q';

  RR = zeros(size(FR,1),size(FR,2));
  GG = zeros(size(FG,1),size(FG,2));
  BB = zeros(size(FB,1),size(FB,2));
% loop over multiple RHS and color components
  for i=1:size(FR,2)
    RR(:,i)=solveRHS(L,U,P,Q1,FR(:,i));  
    GG(:,i)=solveRHS(L,U,P,Q1,FG(:,i));  
    BB(:,i)=solveRHS(L,U,P,Q1,FB(:,i));  
  end
  
end


% Solves single RHS problem for predone LU factorization
function res=solveRHS(L,U,P,Q1,b)
  y1 = L\(P*b);
  y2=U\y1; 
  res=Q1\y2;
end

% Scales [0-1] back to pixel coordinates
function resx=xx(x)
  resx = floor((ix-1)*x+1);
end

% Scales [0-1] back to pixel coordinates
function resy=yy(y)
  resy = floor((iy-1)*y+1);
end


% Helper subroutine for integration over bitmap
function val=bitmp(M,x,y)
  val = zeros(size(x));
  for i=1:size(x,1)
    for j=1:size(x,1)
      val(i,j)=M(xx(x(1,i)),yy(y(1,j)));
    end
  end
end


% Subroutine computing order of polynomials
function p=compute_p(knot_vector)

% first entry in knot_vector
  initial = knot_vector(1);
% lenght of knot_vector
  kvsize = size(knot_vector,2);
  p = 0;

% checking number of repetitions of first entry in knot_vector
  while (p+2 <= kvsize) && (initial == knot_vector(p+2))
    p = p+1;
  end
  
  return
end
  
  
  
  
  
% Subroutine checking sanity of knot_vector
function t=check_sanity(knot_vector,p)

  initial = knot_vector(1);
  kvsize = size(knot_vector,2);

  t = true;
  counter = 1;

% if number of repeated knots at the beginning of knot_vector doesn't match polynomial order
  for i=1:p+1
    if (initial ~= knot_vector(i))
% return FALSE
      t = false;
      return
    end
  end

% if there are too many repeated knots in the middle of knot_vector
  for i=p+2:kvsize-p-1
    if (initial == knot_vector(i))
      counter = counter + 1;
      if (counter > p)
% return FALSE
        t = false;
        return
      end
    else
      initial = knot_vector(i);
      counter = 1;
    end
  end

  initial = knot_vector(kvsize);
  
% if number of repeated knots at the end of knot_vector doesn't match polynomial order
  for i=kvsize-p:kvsize
    if (initial ~= knot_vector(i))
% return FALSE
      t = false;
      return
    end
  end
  
% if subsequent element in knot_vector is smaller than previous one
  for i=1:kvsize-1
    if (knot_vector(i)>knot_vector(i+1))
% return FALSE
      t = false;
      return
    end
  end

  return
end


% Subroutine computing basis functions according to recursive Cox-de-Boor formulae
function y=compute_spline(knot_vector,p,nr,x)
  
% function (x-x_i)/(x_{i-p}-x_i)
  fC= @(x,a,b) (x-a)/(b-a);
% function (x_{i+p+1}-x)/(x_{i+p+1}-x_{i+1})
  fD= @(x,c,d) (d-x)/(d-c);
  
% x_i  
  a = knot_vector(nr);
% x_{i-p} 
  b = knot_vector(nr+p);
% x_{i+1}
  c = knot_vector(nr+1);
% x_{i+p+1}
  d = knot_vector(nr+p+1);

% linear function for p=0
  if (p==0)
    y = 0 .* (x < a) + 1 .* (a <= x & x <= d) + 0 .* (x > d);
    return
  end

% B_{i,p-1}  
  lp = compute_spline(knot_vector,p-1,nr,x);
% B_{i+1,p-1}
  rp = compute_spline(knot_vector,p-1,nr+1,x);
  
% (x-x_i)/(x_{i-p)-x_i)*B_{i,p-1}  
  if (a==b)
% if knots in knot_vector are repeated we have to include it in formula
    y1 = 0 .* (x < a) + 1 .* (a <= x & x <= b) + 0 .* (x > b);
  else
    y1 = 0 .* (x < a) + fC(x,a,b) .* (a <= x & x <= b) + 0 .* (x > b);
  end
  
% (x_{i+p+1}-x)/(x_{i+p+1)-x_{i+1})*B_{i+1,p-1}
  if (c==d)
% if knots in knot_vector are repeated we have to include it in formula
    y2 = 0 .* (x < c) + 1 .* (c < x & x <= d) + 0 .* (d < x);
  else
    y2 = 0 .* (x < c) + fD(x,c,d) .* (c < x & x <= d) + 0 .* (d < x);
  end
  
  y = lp .* y1 + rp .* y2;
  return
end

  
  
% Computes number of elements in given knot_vector
function n=number_of_elements(knot_vector,p)
  initial = knot_vector(1);
  kvsize = size(knot_vector,2);
  n = 0;
  
  for i=1:kvsize-1
    if (knot_vector(i) ~= initial)
      initial = knot_vector(i);
      n = n+1;
    end
  end
end
  
% Creates simple knot_vector without repetitions in the middle
function knot=simple_knot(elems, p)
  pad = ones(1, p);
  knot = [0 * pad, 0:elems, elems * pad];
  knot = knot/elems;
end


% Computes number of degrees of freedom over give knot_vector
function n=number_of_dofs(knot,p)
  n = length(knot) - p - 1;
end


% Finds index of first knot in knot_vector related to give element
function first=first_dof_on_element(knot_vector,p,elem_number)
 [l,h] = element_boundary(knot_vector,p,elem_number);
 first = find(knot_vector==l, 1, 'last') - p;
end
  

% Finds lower and higher boundary of element
function [low,high]=element_boundary(knot_vector,p,elem_number)
  initial = knot_vector(1);
  kvsize = size(knot_vector,2);
  k = 0;
  low=0;
  high=0;
  
  for i=1:kvsize
    if (knot_vector(i) ~= initial)
      initial = knot_vector(i);
      k = k+1;
    end
    if (k == elem_number)
      low = knot_vector(i-1);
      high = knot_vector(i);
      return;
    end
  end
end


% Returns range (indexes) of nonzero functions over element on given knot_vector
function [low,high]=dofs_on_element(knot_vector,p,elem_number)
  low = first_dof_on_element(knot_vector,p,elem_number);
% we expect exactly p+1 nonzero functions over element
  high = low + p;
end


% Row vector of points of the k-point Gaussian quadrature on [a, b]
function xs=quad_points(a, b, k)
% mapping points
  map = @(x) 0.5 * (a * (1 - x) + b * (x + 1));
  
  switch (k)
    case 1
      xs = [0];
    case 2
      xs = [-1/sqrt(3), ...
             1/sqrt(3)];
    case 3
      xs = [-sqrt(3/5), ...
             0,         ...
             sqrt(3/5)];
    case 4
      xs = [-sqrt((3+2*sqrt(6/5))/7), ...
             sqrt((3-2*sqrt(6/5))/7), ...
             sqrt((3-2*sqrt(6/5))/7), ...
             sqrt((3+2*sqrt(6/5))/7)];
    case 5
      xs = [-1/3*sqrt(5+2*sqrt(10/7)), ...
            -1/3*sqrt(5-2*sqrt(10/7)), ...
             0,                        ...
             1/3*sqrt(5-2*sqrt(10/7)), ...
             1/3*sqrt(5+2*sqrt(10/7))];
    otherwise
      xs = [-1/3*sqrt(5+2*sqrt(10/7)), ...
            -1/3*sqrt(5-2*sqrt(10/7)), ...
             0,                        ...
             1/3*sqrt(5-2*sqrt(10/7)), ...
             1/3*sqrt(5+2*sqrt(10/7))];
  end
  xs = map(xs);
end

% Row vector of weights of the k-point Gaussian quadrature on [a, b]
function ws=quad_weights(a, b, k)
  switch (k)
    case 1
      ws = [2];
    case 2
      ws = [1, 1];
    case 3
      ws = [5/9, ...
            8/9, ...
            5/9];
    case 4
      ws = [(18-sqrt(30))/36, ...
            (18+sqrt(30))/36, ...
            (18+sqrt(30))/36, ...
            (18-sqrt(30))/36];
    case 5
      ws = [(322-13.0*sqrt(70))/900, ...
            (322+13.0*sqrt(70))/900, ...
            128/225,                 ...
            (322+13.0*sqrt(70))/900, ...
            (322-13.0*sqrt(70))/900];
    otherwise
      ws = [(322-13.0*sqrt(70))/900, ...
            (322+13.0*sqrt(70))/900, ...
            128/225,                 ...
            (322+13.0*sqrt(70))/900, ...
            (322-13.0*sqrt(70))/900];
  end
end

end