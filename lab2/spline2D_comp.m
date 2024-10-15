function spline2D_comp()
spline2Duniform()
%spline2Dmatrix()
return

function t=check_sanity(knot_vector,p)

  initial = knot_vector(1);
  kvsize = size(knot_vector,2);

  t = true;
  counter = 1;

  for i=1:p+1
    if (initial ~= knot_vector(i))
      t = false;
      return
    end
  end


  for i=p+2:kvsize-p-1
    if (initial == knot_vector(i))
      counter = counter + 1;
      if (counter > p)
        t = false  ;
        return
      end
    else
      initial = knot_vector(i);
      counter = 1;
    end
  end

  initial = knot_vector(kvsize);

  for i=kvsize-p:kvsize
    if (initial ~= knot_vector(i))
      t = false;
      return
    end
  end
  
  for i=1:kvsize-1
    if (knot_vector(i)>knot_vector(i+1))
      t = false;
    end
  end

  return

end


function spline2Dmatrix()
knot_vectorx=[0 1 2];

knot_vectory=[0 1 2];

precision = 0.1;

%macros
compute_nr_basis_functions = @(knot_vector,p) size(knot_vector, 2) - p - 1;
mesh   = @(a,c) [a:precision*(c-a):c];

%splines in x
px = compute_p(knot_vectorx);
tx = check_sanity(knot_vectorx,px);
nrx = compute_nr_basis_functions(knot_vectorx,px);

x_begin = knot_vectorx(1);
x_end = knot_vectorx(size(knot_vectorx,2));

x=mesh(x_begin,x_end);

%splines in y
py = compute_p(knot_vectory);
ty = check_sanity(knot_vectory,py);
nry = compute_nr_basis_functions(knot_vectorx,py);

y_begin = knot_vectory(1);
y_end = knot_vectory(size(knot_vectory,2));

y=mesh(y_begin,y_end);

%X and Y coordinates of points over the 2D mesh
[X,Y]=meshgrid(x,y);

M(1:nrx*nry,1:nrx*nry)=0.0;
hold on
for i=1:nrx
  %compute values of 
  vx=compute_spline(knot_vectorx,px,i,X);
  for j=1:nry
    vy=compute_spline(knot_vectory,py,j,Y);
    spline1=vx.*vy;
    for k=1:nrx
      %compute values of 
      vx2=compute_spline(knot_vectorx,px,k,X);
      for l=1:nry
        vy2=compute_spline(knot_vectory,py,l,Y);
        spline2=vx2.*vy2;
        crossection = spline1.*spline2;
        M(i+(j-1)*nry,k+(l-1)*nry)=M(i+(j-1)*nry,k+(l-1)*nry)+norm(crossection);
      end
    end
  end
end
crossection;
hold off

end

function spline2Duniform()
knot_vectorx=[0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 29];
knot_vectory=[0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 29];
coeff_vector=[
[1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1];
[0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1];
[1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1];
[1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1];
[1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0];
[0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0];
[1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0];
[1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0];
[0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1];
[0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0];
[1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1];
[1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1];
[1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1];
[0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1];
[1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1];
[0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1];
[0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0];
[0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1];
[1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1];
[0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0];
[1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1];
[1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1];
[1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1];
[1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1];
[0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1];
[1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0];
[1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1];
[1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1];
[1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1];
[1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1]
];

precision = 0.01;

%macros
compute_nr_basis_functions = @(knot_vector,p) size(knot_vector, 2) - p - 1
mesh   = @(a,c) [a:precision*(c-a):c];

%splines in x
px = compute_p(knot_vectorx);
tx = check_sanity(knot_vectorx,px);
nrx = compute_nr_basis_functions(knot_vectorx,px);

x_begin = knot_vectorx(1);
x_end = knot_vectorx(size(knot_vectorx,2));

x=mesh(x_begin,x_end);

%splines in y
py = compute_p(knot_vectorx);
ty = check_sanity(knot_vectorx,py);
nry = compute_nr_basis_functions(knot_vectorx,py);

y_begin = knot_vectory(1);
y_end = knot_vectory(size(knot_vectory,2));

y=mesh(y_begin,y_end);

%X and Y coordinates of points over the 2D mesh
[X,Y]=meshgrid(x,y);

% BEGIN changes
u = 0;

hold on
for i=1:nrx
  %compute values of 
  vx=compute_spline(knot_vectorx,px,i,X);
  for j=1:nry
    vy=compute_spline(knot_vectory,py,j,Y);
    u = u + vx.*vy.*coeff_vector(i,j);
  end
end
surf(X,Y,u);
hold off
% END changes

end

function spline()

knot_vector=[0 0 0  1 1 2 2 2 ];

precision = 0.0001;


compute_nr_basis_functions = @(knot_vector,p) size(knot_vector, 2) - p - 1;
mesh   = @(a,c) [a:precision*(c-a):c];


p = compute_p(knot_vector);
t = check_sanity(knot_vector,p);
nr = compute_nr_basis_functions(knot_vector,p);

x_begin = knot_vector(1);
x_end = knot_vector(size(knot_vector,2));


x=mesh(x_begin,x_end);


y=compute_spline(knot_vector,p,1,x);
plot(x,y)

hold on

for i=2:nr
  y=compute_spline(knot_vector,p,i,x);
  plot(x,y)
end

hold off


end


function p=compute_p(knot_vector)

  initial = knot_vector(1);
  kvsize = size(knot_vector,2);
  i=1;

  while (i+1 < kvsize) && (initial == knot_vector(i+1))
    i=i+1;
  end
  
  p = i-1;
  
  return
end



function y=compute_spline(knot_vector,p,nr,x)
  
  fC= @(x,a,b) (x)/(b-a)-a/(b-a);
  fD= @(x,c,d) (1-x)/(d-c)+(d-1)/(d-c);
  
  
  a = knot_vector(nr);
  b = knot_vector(nr+p);
  c = knot_vector(nr+1);
  d = knot_vector(nr+p+1);

  if (p==0)
    y = 0 .* (x < a) + 1 .* (a <= x & x <= d) + 0 .* (x > d);
    return
  end
  
  lp = compute_spline(knot_vector,p-1,nr,x);
  rp = compute_spline(knot_vector,p-1,nr+1,x);
  
  if (a==b)
    y1 = 0 .* (x < a) + 1 .* (a <= x & x <= b) + 0 .* (x > b);
  else
    y1 = 0 .* (x < a) + fC(x,a,b) .* (a <= x & x <= b) + 0 .* (x > b);
  end
  
  if (c==d)
    y2 = 0 .* (x < c) + 1 .* (c < x & x <= d) + 0 .* (d < x);
  else
    y2 = 0 .* (x < c) + fD(x,c,d) .* (c < x & x <= d) + 0 .* (d < x);
  end
  
  y = lp .* y1 + rp .* y2;
  return
  
end

end