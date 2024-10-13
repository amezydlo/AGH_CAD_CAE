% This subroutine draws B-Spline basis functions using given knot vector.
% How to run
%
% splines(number of points to draw, knot vector as a matrix)
%
% Examples
%
% draws evenly distributed knots, 1-st order splines, 3 elements
% splines (1000, [1,1,2,3,4,4])
%
% draws unevenly distributed knots, 3-rd order splines, 6 elements
% splines (100000,[0,0,0,0,1,3,6,10,11,11.5,11.5,11.5,11.5])
%
% draws evenly distributed knots with repeated knots, 2-nd order splines, 3 elements
% splines (1000, [1,1,1,2,2,3,4,4,4])

% splines_comb(1000, [0,0,0,1,1,2,3,4,5,6,7,8,8,8], [1,1,1,1,1,1,1,1,1,1,1]);

% p=2 version
% splines_comb(1000, [0, 0, 0, 0.62, 3.44, 4.38, 4.84, 5.16, 5.31, 7.81, 7.97, 9.38, 10.0, 10.78, 11.56, 13.12, 13.44, 15.0, 15.62, 16.09, 17.34, 18.12, 19.22, 19.53, 20.62, 21.41, 21.72, 22.5, 22.81, 25.0, 25.62, 25.94, 26.09, 27.97, 28.28, 28.59, 29.22, 30.0, 31.87, 32.03, 33.28, 33.75, 34.69, 34.84, 36.25, 37.03, 37.19, 38.59, 39.53, 41.25, 42.97, 43.91, 45.0, 45.16, 46.25, 48.44, 48.75, 49.06, 51.56, 52.03, 52.34, 52.81, 53.59, 55.0, 55.16, 56.41, 57.19, 58.44, 60.62, 61.72, 62.97, 63.12, 65.78, 65.94, 67.34, 69.06, 70.16, 72.19, 72.97, 73.12, 73.59, 74.38, 75.16, 75.31, 75.94, 76.25, 76.72, 77.34, 77.66, 79.06, 79.53, 81.56, 82.19, 82.97, 84.53, 86.41, 88.28, 90.78, 92.03, 92.81, 93.59, 94.22, 95.0, 96.25, 97.81, 98.12, 99.38, 100, 100, 100], [39.22, 39.22, 39.38, 39.22, 35.0, 34.84, 37.5, 37.5, 36.41, 36.25, 38.28, 38.28, 38.75, 38.75, 37.81, 37.66, 38.59, 39.38, 39.84, 39.38, 39.06, 38.59, 38.44, 37.97, 36.09, 36.09, 36.56, 36.72, 36.72, 37.81, 42.66, 42.66, 42.03, 40.62, 40.47, 40.16, 40.0, 39.06, 38.91, 35.62, 35.62, 37.97, 38.12, 38.44, 39.06, 39.06, 39.84, 40.16, 40.47, 38.12, 36.72, 38.12, 39.22, 37.97, 38.59, 44.06, 43.75, 40.0, 37.03, 41.72, 42.66, 41.87, 39.53, 39.53, 38.91, 38.75, 38.91, 36.56, 36.88, 37.81, 37.81, 36.72, 36.25, 36.25, 36.25, 35.47, 35.78, 36.56, 36.09, 35.78, 37.03, 37.03, 36.09, 35.94, 36.56, 37.97, 37.03, 35.62, 35.47, 36.72, 37.5, 35.94, 35.62, 36.41, 36.72, 36.72, 36.88, 37.34, 37.66, 37.97, 37.66, 36.88, 37.66, 38.44, 38.59, 38.28, 38.29])

% p=1 version - better
% splines_comb(1000, [0, 0, 0.62, 3.44, 4.38, 4.84, 5.16, 5.31, 7.81, 7.97, 9.38, 10.0, 10.78, 11.56, 13.12, 13.44, 15.0, 15.62, 16.09, 17.34, 18.12, 19.22, 19.53, 20.62, 21.41, 21.72, 22.5, 22.81, 25.0, 25.62, 25.94, 26.09, 27.97, 28.28, 28.59, 29.22, 30.0, 31.87, 32.03, 33.28, 33.75, 34.69, 34.84, 36.25, 37.03, 37.19, 38.59, 39.53, 41.25, 42.97, 43.91, 45.0, 45.16, 46.25, 48.44, 48.75, 49.06, 51.56, 52.03, 52.34, 52.81, 53.59, 55.0, 55.16, 56.41, 57.19, 58.44, 60.62, 61.72, 62.97, 63.12, 65.78, 65.94, 67.34, 69.06, 70.16, 72.19, 72.97, 73.12, 73.59, 74.38, 75.16, 75.31, 75.94, 76.25, 76.72, 77.34, 77.66, 79.06, 79.53, 81.56, 82.19, 82.97, 84.53, 86.41, 88.28, 90.78, 92.03, 92.81, 93.59, 94.22, 95.0, 96.25, 97.81, 98.12, 99.38, 100, 100], [39.22, 39.22, 39.38, 39.22, 35.0, 34.84, 37.5, 37.5, 36.41, 36.25, 38.28, 38.28, 38.75, 38.75, 37.81, 37.66, 38.59, 39.38, 39.84, 39.38, 39.06, 38.59, 38.44, 37.97, 36.09, 36.09, 36.56, 36.72, 36.72, 37.81, 42.66, 42.66, 42.03, 40.62, 40.47, 40.16, 40.0, 39.06, 38.91, 35.62, 35.62, 37.97, 38.12, 38.44, 39.06, 39.06, 39.84, 40.16, 40.47, 38.12, 36.72, 38.12, 39.22, 37.97, 38.59, 44.06, 43.75, 40.0, 37.03, 41.72, 42.66, 41.87, 39.53, 39.53, 38.91, 38.75, 38.91, 36.56, 36.88, 37.81, 37.81, 36.72, 36.25, 36.25, 36.25, 35.47, 35.78, 36.56, 36.09, 35.78, 37.03, 37.03, 36.09, 35.94, 36.56, 37.97, 37.03, 35.62, 35.47, 36.72, 37.5, 35.94, 35.62, 36.41, 36.72, 36.72, 36.88, 37.34, 37.66, 37.97, 37.66, 36.88, 37.66, 38.44, 38.59, 38.28])
     

function splines_comp(precision, knot_vector, coeff_vector)

% subroutine calculating number of basis functions
compute_nr_basis_functions = @(knot_vector,p) size(knot_vector, 2) - p - 1;
% subroutine generating mesh for drawing basis functions
mesh = @(a,c) [a:(c-a)/precision:c];
% subroutine drawing basis functions
% plot_spline = @(knot_vector,p,nr,x) plot(x,compute_splines(knot_vector,p,nr,x));

% computing order of polynomials
p = compute_p(knot_vector);
% validation of knot vector construction
t = check_sanity(knot_vector,p);
% if knot vector is poorly constructed - stop further processing
if (~t)
  disp("Poorly constructed knot_vector")
  return
end
% computating number of basis functions
nr = compute_nr_basis_functions(knot_vector,p);

% beginning of drawing range
x_begin = knot_vector(1);
% end of drawing range
x_end = knot_vector(size(knot_vector,2));
x=mesh(x_begin,x_end);


% My implementation
if (nr ~= size(coeff_vector, 2)) 
    disp("Poorly constructed coeff_vector")
end

u = zeros(size(x));
for j=1:nr
    u = u + coeff_vector(j) * compute_splines(knot_vector, p, j, x);
end

% Plot with background image
figure;
plot(x, u, 'LineWidth', 3, 'color', 'red');
axis tight;
hold on;
I = imread('image.jpg'); % Load image
I_flipped = flipud(I);    % Flip image vertically - to fix coordinates
h = image('CData', I_flipped, 'XData', [0, 100], 'YData', [0, 100]);
uistack(h, 'bottom');     % Send image to background
hold off;
saveas(gcf, 'spline_on_image.png');

% Plot without background image
figure;
plot(x, u, 'LineWidth', 3, 'color', 'red');
axis tight;
xlim([0, 100]);
ylim([0, 100]);  
saveas(gcf, 'spline_only.png');

% End of my implementation



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

% Subroutine computing basis functions according to recursive Cox-de-Boor formulae
function y=compute_splines(knot_vector,p,nr,x)

y=compute_spline(knot_vector,p,nr,x);
%y=y.*sin(pi/10*x);
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



end