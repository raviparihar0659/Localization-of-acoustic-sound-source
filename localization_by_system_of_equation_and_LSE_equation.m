global c s1 s2 s3 s4 s5 s6 s7 s8 Ts;




%%

s1=[0,0,0];s2=[1,0,0];s0=[30,20,25];
% tau= abs(norm(s0-s1)-norm(s0-s2))/343;

% Define the equation as a function handle
f = @(x, y, z) sqrt((x-1).^2 + (y).^2 + (z).^2) - sqrt(x.^2 + y.^2 + z.^2) - (norm(s0-s2)-norm(s0-s1));

% Create a finer grid of points
[x, y, z] = meshgrid(-35:0.5:35, -35:0.5:35, -35:0.5:35);

% Evaluate the function on the grid
v = f(x, y, z);

% Create an isosurface for the value 0
isosurface_handle = isosurface(x, y, z, v, 0);

% Check if isosurface is empty
if ~isempty(isosurface_handle.vertices)
    % Plot the isosurface
    p = patch(isosurface_handle);
    isonormals(x, y, z, v, p)
    set(p, 'FaceColor', 'red', 'EdgeColor', 'none');

    % Add lighting
    camlight; 
    lighting phong;

    % Label the axes
    xlabel('x');
    ylabel('y');
    zlabel('z');

    % Set the view and add grid
    view(3);
    grid on;

    % Add title
    title('Plot of \sqrt{(x-1)^2 + (y-4)^2 + (z-2)^2} - \sqrt{x^2 + y^2 + z^2} = 5');
else
    disp('No isosurface found within the grid range.');
end


%%

s1=[0,0,0];s2=[1.1,0,0];s3=[1,0.1,0];s4=[1,0,0.1];  s0=[30,20,25];
% Define the equations as function handles
f1 = @(x, y, z) sqrt((x-1.1).^2 + (y).^2 + (z).^2) - sqrt(x.^2 + y.^2 + z.^2) - (norm(s0-s2)-norm(s0-s1));
f2 = @(x, y, z) sqrt((x-1).^2 + (y-0.1).^2 + (z).^2) - sqrt(x.^2 + y.^2 + z.^2) - (norm(s0-s3)-norm(s0-s1));
f3 = @(x, y, z) sqrt((x-1).^2 + (y).^2 + (z-0.1).^2) - sqrt(x.^2 + y.^2 + z.^2) - (norm(s0-s4)-norm(s0-s1));


% Create a finer grid of points
[x, y, z] = meshgrid(-50:0.5:50, -50:0.5:50, -50:0.5:50);

% Evaluate the functions on the grid
v1 = f1(x, y, z);
v2 = f2(x, y, z);
v3 = f3(x, y, z);

% Create figure
figure;

% Plot the first isosurface
isosurface_handle1 = isosurface(x, y, z, v1, 0);
p1 = patch(isosurface_handle1);
isonormals(x, y, z, v1, p1)
set(p1, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

hold on;

% Plot the second isosurface
isosurface_handle2 = isosurface(x, y, z, v2, 0);
p2 = patch(isosurface_handle2);
isonormals(x, y, z, v2, p2)
set(p2, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Plot the third isosurface
isosurface_handle3 = isosurface(x, y, z, v3, 0);
p3 = patch(isosurface_handle3);
isonormals(x, y, z, v3, p3)
set(p3, 'FaceColor', 'green', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Add lighting
camlight; 
lighting phong;

% Label the axes
xlabel('x');
ylabel('y');
zlabel('z');

% Set the view and add grid
view(3);
grid on;

% Add title
title('Plot of Multiple Implicit Surfaces');
hold off;

%%

s1=[0,0,0];s2=[1.1,0,0];s3=[1,0.1,0];s4=[1,0,0.1];s5=[1,0,0];  s0=[30,20,50];
system_of_equations = @(vars) [
    sqrt((vars(1)-1.1).^2 + (vars(2)).^2 + (vars(3)).^2) - sqrt(vars(1).^2 + vars(2).^2 + vars(3).^2) - (norm(s0-s2)-norm(s0-s1));
    sqrt((vars(1)-1).^2 + (vars(2)-0.1).^2 + (vars(3)).^2) - sqrt(vars(1).^2 + vars(2).^2 + vars(3).^2) - (norm(s0-s3)-norm(s0-s1));
    sqrt((vars(1)-1).^2 + (vars(2)).^2 + (vars(3)-0.1).^2) - sqrt(vars(1).^2 + vars(2).^2 + vars(3).^2) - (norm(s0-s4)-norm(s0-s1));
    sqrt((vars(1)-1).^2 + (vars(2)).^2 + (vars(3)).^2) - sqrt(vars(1).^2 + vars(2).^2 + vars(3).^2) - (norm(s0-s5)-norm(s0-s1));
];
% Initial guess for the solver
initial_guess = [0, 0, 0];
% Solve the system of equations
options = optimoptions('fsolve', 'Display', 'iter');
solution = fsolve(system_of_equations, initial_guess, options);
% Display the solution
disp('Solution:');
disp(solution);
%%

s1=[0,0,0];s2=[0.5,0,0];s3=[0,0.5,0];s4=[0,0,0.5];s5=[-0.5,0,0]; 
% here s1 is reference
s0=[60,45,50]; % original source location
c=343;
fs=64000*32; % for correct solution it requires very high sampling frequency ( need to upsample by factor 32 is sampling frequency is 64000Hz)
system_of_equations = @(vars) [
    sqrt((vars(1)-s2(1)).^2 + (vars(2)-s2(2)).^2 + (vars(3)-s2(3)).^2) - sqrt(vars(1).^2 + vars(2).^2 + vars(3).^2) - c*(round(((norm(s0-s2)-norm(s0-s1))/c)*fs)/fs);
    sqrt((vars(1)-s3(1)).^2 + (vars(2)-s3(2)).^2 + (vars(3)-s3(3)).^2) - sqrt(vars(1).^2 + vars(2).^2 + vars(3).^2) - c*(round(((norm(s0-s3)-norm(s0-s1))/c)*fs)/fs);
    sqrt((vars(1)-s4(1)).^2 + (vars(2)-s4(2)).^2 + (vars(3)-s4(3)).^2) - sqrt(vars(1).^2 + vars(2).^2 + vars(3).^2) - c*(round(((norm(s0-s4)-norm(s0-s1))/c)*fs)/fs);
    sqrt((vars(1)-s5(1)).^2 + (vars(2)-s5(2)).^2 + (vars(3)-s5(3)).^2) - sqrt(vars(1).^2 + vars(2).^2 + vars(3).^2) - c*(round(((norm(s0-s5)-norm(s0-s1))/c)*fs)/fs);
];

% Initial guess for the solver
initial_guess = [0, 0, 0];

% Solve the system of equations
options = optimoptions('fsolve', 'Display', 'iter');
solution = fsolve(system_of_equations, initial_guess, options);

% Display the solution
disp('Solution:');
disp(solution);

%%
s1 = [0, 0, 0];
s2 = [0.2, 0, 0];
s3 = [0, 0.2, 0];
s4 = [0, 0, 0.2];
sref= [-1,0,0];
s =[43,90,65];
ts= 1/(64000*64);
c=343;

% TDOA=[
%     (norm(s-s1)-norm(s-sref))/c, ...
%     (norm(s-s2)-norm(s-sref))/c,...
%     (norm(s-s3)-norm(s-sref))/c,...
%     (norm(s-s4)-norm(s-sref))/c];
 TDOA= [round(((norm(s-s1)-norm(s-sref))/c)/ts)*ts,round(((norm(s-s2)-norm(s-sref))/c)/ts)*ts,round(((norm(s-s3)-norm(s-sref))/c)/ts)*ts,round(((norm(s-s4)-norm(s-sref))/c)/ts)*ts];
% fprintf('%X.Yf', TDOA);
% disp(TDOA);
A= 2*[
    s1-sref ,c*TDOA(1);
    s2-sref , c*TDOA(2);
    s3-sref , c*TDOA(3);
    s4-sref , c*TDOA(4)
    ];
B= [
  norm(s1)^2-norm(sref)^2- (c*TDOA(1))^2;
  norm(s2)^2-norm(sref)^2- (c*TDOA(2))^2;
  norm(s3)^2-norm(sref)^2- (c*TDOA(3))^2;
  norm(s4)^2-norm(sref)^2- (c*TDOA(4))^2
];

XLos = (A' * A) \ (A' * B);

drone_position = XLos(1:3);
fprintf('Estimated Drone Position: X = %.2f, Y = %.2f, Z = %.2f\n', drone_position);



%% writing error equation of localization and searching least square error equation soltion in 3D grid

% now we come up new method  other than all the solutions we solved yet. we
% take 4 microphones , using the least squared error method we solve the
% equations of their path , time difference of arrival will be their cross
% coupled pairs, there is no perticular reference here


% Sensor locations
s1 = [0, 0, 0];
s2 = [1, 0, 0];
s3 = [0, 1, 0];
s4 = [0, 0, 1];

% Assumed source location
s0 = [50, 43, 81];

c=343;

% Define the squared error function as a function handle
f = @(x, y, z) (sqrt((x-1).^2 + y.^2 + z.^2) - sqrt(x.^2 + y.^2 + z.^2) - (norm(s0 - s2) - norm(s0 - s1)) ).^2 + ...
               (sqrt(x.^2 + (y-1).^2 + z.^2) - sqrt(x.^2 + y.^2 + z.^2) - (norm(s0 - s3) - norm(s0 - s1))).^2 + ...
               (sqrt(x.^2 + y.^2 + (z-1).^2) - sqrt(x.^2 + y.^2 + z.^2) - (norm(s0 - s4) - norm(s0 - s1))).^2 + ...
               (sqrt(x.^2 + (y-1).^2 + z.^2) - sqrt((x-1).^2 + y.^2 + z.^2) - (norm(s0 - s3) - norm(s0 - s2))).^2 + ...
               (sqrt(x.^2 + y.^2 + (z-1).^2) - sqrt((x-1).^2 + y.^2 + z.^2) - (norm(s0 - s4) - norm(s0 - s2))).^2 + ...
               (sqrt(x.^2 + y.^2 + (z-1).^2) - sqrt(x.^2 + (y-1).^2 + z.^2) - (norm(s0 - s4) - norm(s0 - s3))).^2;

% Define the grid range and step size
[x, y, z] = meshgrid(45:0.1:55, 40:0.1:50, 78:0.1:85);

% Evaluate the function on the grid
v = f(x, y, z);

% Find the minimum value of v and its corresponding coordinates
[min_v, min_idx] = min(v(:));
[min_x, min_y, min_z] = ind2sub(size(v), min_idx);

% Display the estimated source location
estimated_source = [x(min_x, min_y, min_z), y(min_x, min_y, min_z), z(min_x, min_y, min_z)];
disp(['Estimated source location: ', num2str(estimated_source)]);



% Create figure
figure;

% Plot the isosurface at the minimum value of v
isosurface_value = min_v;
isosurface_handle = isosurface(x, y, z, v, isosurface_value);

% Check if the isosurface handle is empty
if isempty(isosurface_handle.vertices)
    warning('No isosurface found for the value %f.', isosurface_value);
else
    p = patch(isosurface_handle);
    isonormals(x, y, z, v, p);
    set(p, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

    hold on;
end

% Mark the estimated source location
plot3(estimated_source(1), estimated_source(2), estimated_source(3), 'bo', 'MarkerSize', 10, 'LineWidth', 2);

% Adjust visualization properties
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting gouraud

% ploting the squared error:

% Define a slice in the x-y plane at a specific z coordinate (you can adjust z_slice as needed)
z_slice = 30;

% Evaluate f over the slice
f_slice = f(x(:,:,z_slice), y(:,:,z_slice), z(:,:,z_slice));

% Plot the function f
figure;
contourf(x(:,:,z_slice), y(:,:,z_slice), f_slice, 50, 'LineWidth', 0.5);
colorbar;
xlabel('x');
ylabel('y');
title('Slice of f(x, y, z) at z = 25');

hold on;

% Plot the estimated source location
plot(estimated_source(1), estimated_source(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
text(estimated_source(1), estimated_source(2), [' Estimated Source: (' num2str(estimated_source(1)) ', ' num2str(estimated_source(2)) ')']);

% Plot the actual source location
plot(s0(1), s0(2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
text(s0(1), s0(2), [' Actual Source: (' num2str(s0(1)) ', ' num2str(s0(2)) ')']);

hold off;

% Adjust axis limits for clarity
axis equal;



%%
s1 = [0, 0, 0];
s2 = [1, 0, 0];
s3 = [0, 1, 0];
s4 = [0, 0, 1];

% Assumed source location
s0 = [35.53, 37.47, 28.19];
ts=1/(64000*8);
c=343;


% Define the squared error function as a function handle
f = @(x, y, z) (sqrt((x-1).^2 + y.^2 + z.^2) - sqrt(x.^2 + y.^2 + z.^2) - c*round(((norm(s0 - s2) - norm(s0 - s1))/c)/ts)*ts ).^2 + ...
               (sqrt(x.^2 + (y-1).^2 + z.^2) - sqrt(x.^2 + y.^2 + z.^2) - c*round(((norm(s0 - s3) - norm(s0 - s1))/c)/ts)*ts).^2 + ...
               (sqrt(x.^2 + y.^2 + (z-1).^2) - sqrt(x.^2 + y.^2 + z.^2) - c*round(((norm(s0 - s4) - norm(s0 - s1))/c)/ts)*ts).^2 + ...
               (sqrt(x.^2 + (y-1).^2 + z.^2) - sqrt((x-1).^2 + y.^2 + z.^2) - c*round(((norm(s0 - s3) - norm(s0 - s2))/c)/ts)*ts).^2 + ...
               (sqrt(x.^2 + y.^2 + (z-1).^2) - sqrt((x-1).^2 + y.^2 + z.^2) - c*round(((norm(s0 - s4) - norm(s0 - s2))/c)/ts)*ts).^2 + ...
               (sqrt(x.^2 + y.^2 + (z-1).^2) - sqrt(x.^2 + (y-1).^2 + z.^2) - c*round(((norm(s0 - s4) - norm(s0 - s3))/c)/ts)*ts).^2;

% Define the grid range and step size
[x, y, z] = meshgrid(25:0.1:40, 25:0.1:40, 25:0.1:40);

% Evaluate the function on the grid
v = f(x, y, z);

% Find the minimum value of v and its corresponding coordinates
[min_v, min_idx] = min(v(:));
[min_x, min_y, min_z] = ind2sub(size(v), min_idx);

% Display the estimated source location
estimated_source = [x(min_x, min_y, min_z), y(min_x, min_y, min_z), z(min_x, min_y, min_z)];
disp(['Estimated source location: ', num2str(estimated_source)]);



% Create figure
figure;

% Plot the isosurface at the minimum value of v
isosurface_value = min_v;
isosurface_handle = isosurface(x, y, z, v, isosurface_value);

% Check if the isosurface handle is empty
if isempty(isosurface_handle.vertices)
    warning('No isosurface found for the value %f.', isosurface_value);
else
    p = patch(isosurface_handle);
    isonormals(x, y, z, v, p);
    set(p, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

    hold on;
end

% Mark the estimated source location
plot3(estimated_source(1), estimated_source(2), estimated_source(3), 'bo', 'MarkerSize', 10, 'LineWidth', 2);

% Adjust visualization properties
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting gouraud

% ploting the squared error:

% Define a slice in the x-y plane at a specific z coordinate (you can adjust z_slice as needed)
z_slice = 30;

% Evaluate f over the slice
f_slice = f(x(:,:,z_slice), y(:,:,z_slice), z(:,:,z_slice));

% Plot the function f
figure;
contourf(x(:,:,z_slice), y(:,:,z_slice), f_slice, 50, 'LineWidth', 0.5);
colorbar;
xlabel('x');
ylabel('y');
title('Slice of f(x, y, z) at z = 25');

hold on;

% Plot the estimated source location
plot(estimated_source(1), estimated_source(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
text(estimated_source(1), estimated_source(2), [' Estimated Source: (' num2str(estimated_source(1)) ', ' num2str(estimated_source(2)) ')']);

% Plot the actual source location
plot(s0(1), s0(2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
text(s0(1), s0(2), [' Actual Source: (' num2str(s0(1)) ', ' num2str(s0(2)) ')']);

hold off;

% Adjust axis limits for clarity
axis equal;


%%
% Define the grid range and step size
global X Y Z;
[X, Y, Z] = meshgrid(5:1:105, 5:1:105, 30:1:40);

% Function to find the estimated source location
function estimated_source_location = find_estimated_location(s0)
    global X Y Z;
    s1 = [0, 0, 0];
    s2 = [1, 0, 0];
    s3 = [0, 1, 0];
    s4 = [0, 0, 1];

    ts = 1 / (64000 * 8);
    c = 343;

    % Define the squared error function as a function handle
    f = @(x, y, z) (sqrt((x-1).^2 + y.^2 + z.^2) - sqrt(x.^2 + y.^2 + z.^2) - c * round(((norm(s0 - s2) - norm(s0 - s1)) / c) / ts) * ts ).^2 + ...
                   (sqrt(x.^2 + (y-1).^2 + z.^2) - sqrt(x.^2 + y.^2 + z.^2) - c * round(((norm(s0 - s3) - norm(s0 - s1)) / c) / ts) * ts).^2 + ...
                   (sqrt(x.^2 + y.^2 + (z-1).^2) - sqrt(x.^2 + y.^2 + z.^2) - c * round(((norm(s0 - s4) - norm(s0 - s1)) / c) / ts) * ts).^2 + ...
                   (sqrt(x.^2 + (y-1).^2 + z.^2) - sqrt((x-1).^2 + y.^2 + z.^2) - c * round(((norm(s0 - s3) - norm(s0 - s2)) / c) / ts) * ts).^2 + ...
                   (sqrt(x.^2 + y.^2 + (z-1).^2) - sqrt((x-1).^2 + y.^2 + z.^2) - c * round(((norm(s0 - s4) - norm(s0 - s2)) / c) / ts) * ts).^2 + ...
                   (sqrt(x.^2 + y.^2 + (z-1).^2) - sqrt(x.^2 + (y-1).^2 + z.^2) - c * round(((norm(s0 - s4) - norm(s0 - s3)) / c) / ts) * ts).^2;

    % Evaluate the function on the grid
    v = f(X, Y, Z);

    % Find the minimum value of v and its corresponding coordinates
    [~, min_idx] = min(v(:));
    [min_x, min_y, min_z] = ind2sub(size(v), min_idx);

    % Display the estimated source location
    estimated_source_location = [X(min_x, min_y, min_z), Y(min_x, min_y, min_z), Z(min_x, min_y, min_z)];
    
end

xcor = linspace(10, 100, 91);
ycor = linspace(10, 100, 91);
z = 35;
array3D = zeros(length(xcor), length(ycor), 1);

for i = 1:length(xcor)
    for j = 1:length(ycor)
        s = [xcor(i), ycor(j), z];
        
        estimated_source_location = find_estimated_location(s);
        % disp(s);
        % disp(S_est');
        array3D(i, j, :) = abs(norm(estimated_source_location - s));
    end
end
%%

% Plot the data
[X_grid, Y_grid] = meshgrid(xcor, ycor);
Z_grid = squeeze(array3D(:, :, 1));

figure;
surf(X_grid, Y_grid, Z_grid);
xlabel('X Position');
ylabel('Y Position');
zlabel('Error');
title('Localization Error in 3D Space');
colorbar;


