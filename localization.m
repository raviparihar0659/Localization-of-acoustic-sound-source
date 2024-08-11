%%  use of optimization tool to solve the squared error; first we assumed that we have exact value of time difference of arrival
%   there is need of optimization toolbox install in matlab
% Define the coordinates of the ith mic location (ai, bi, ci)


% S is matrix formed from loactions of microphones
S = [ 
     1, 0, 0;
     0, 1, 0;
     0, 0, 1;
     1, 1, 1
];
s =[58.7,35,83.5]; %source location

% Define vector d
d = [
    norm(s-S(1,:))-norm(s);
    norm(s-S(2,:))-norm(s);
    norm(s-S(3,:))-norm(s);
    norm(s-S(4,:))-norm(s)
];

% Calculate vector c (squared norms of the rows of S)
c = sum(S.^2, 2); % N x 1 vector

% Define the objective function to minimize
objective = @(y) sum((S * y + d * norm(y) - 0.5 * (c - d.^2)).^2);

% Initial guess for y
y0 = rand(size(S, 2), 1);

% Use MATLAB's optimization function to find y that minimizes the objective function
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter');
[y_opt, fval] = fminunc(objective, y0, options);

% Display the results
disp('Optimal y:');
disp(y_opt);
disp('Minimum squared error:');
disp(fval);

%if we have exact value of TDOA then we can esitimate the location with
%error of within 0.02 percent
% for the original location of source s =[58.7,35,83.5] , estimated value
% of source is found
%Optimal y:
   % 58.6892
   % 34.9936
   % 83.4846

%% use of optimization tool to solve the squared error; now we assumed that we have quentized value of time difference of arrival(estimated TDOA will be integer multiple of sampling time)
% Define the coordinates of the ith microphone location (ai, bi, ci)


S = [ 
     1, 0, 0;
     0, 1, 0;
     0, 0, 1;
     1, 1, 1
];

s =[80,35.9,70.5];

ts= 1/(64000*8);
v=343;
% Define vector d
d = [
    v*round(((norm(s-S(1,:))-norm(s))/v)/ts)*ts;
    v*round(((norm(s-S(2,:))-norm(s))/v)/ts)*ts;
    v*round(((norm(s-S(3,:))-norm(s))/v)/ts)*ts;
    v*round(((norm(s-S(4,:))-norm(s))/v)/ts)*ts
];

% Calculate vector c (squared norms of the rows of S)
c = sum(S.^2, 2); % N x 1 vector

% Define the objective function to minimize
objective = @(y) sum((S * y + d * norm(y) - 0.5 * (c - d.^2)).^2);

% Initial guess for y
y0 = rand(size(S, 2), 1);

% Use MATLAB's optimization function to find y that minimizes the objective function
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter');
[y_opt, fval] = fminunc(objective, y0, options);

% Display the results
disp('Optimal y:');
disp(y_opt);
disp('Minimum squared error:');
disp(fval);

%%
% Microphone Array Configuration:
% The matrix S defines the positions of six microphones arranged in a symmetric pattern along the axes. 
% Each row represents the coordinates of a microphone in 3D space.

%TDOA Calculation:
% For a given source position s, the code calculates the TDOA for each microphone relative to an assumed reference point.
% The TDOA is quantized using a sampling period ts, which is derived from the specified sampling frequency.

% a microphone s0 is located at origin which is worked as reference
% microphone

% when we ant to calculate the localization of a sound source; we need to calculate first tdoa t10,t20,t30,t40,t50,t60;
% for given microphones location 
% matrix S = [s1;
%             s2;
%             s3;
%             s4;
%             s5;
%             s6] 
% matrix d = [c*t10;
%             c* t20;
%             c*t30;
%             c*t40;
%             c*t50;
%             c*t60]
% from these S and d matrices apply as give in below

function locate =find(s)

    S = [ 
         1, 0, 0;
         0, 1, 0;
         0, 0, 1;
         -1, 0, 0;
         0, -1, 0;
         0, 0, -1   
    ];

    ts= 1/(64000*8*4);
    v=343;

    d = [
        v*round(((norm(s-S(1,:))-norm(s))/v)/ts)*ts;
        v*round(((norm(s-S(2,:))-norm(s))/v)/ts)*ts;
        v*round(((norm(s-S(3,:))-norm(s))/v)/ts)*ts;
        v*round(((norm(s-S(4,:))-norm(s))/v)/ts)*ts;
        v*round(((norm(s-S(5,:))-norm(s))/v)/ts)*ts;
        v*round(((norm(s-S(6,:))-norm(s))/v)/ts)*ts    
    ];
    
    c = sum(S.^2, 2); % N x 1 vector

    objective = @(y) sum((S * y + d * norm(y) - 0.5 * (c - d.^2)).^2);

    y0 = rand(size(S, 2), 1); 

    % Use MATLAB's optimization function to find y that minimizes the objective function
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter');
    [y_opt, ~] = fminunc(objective, y0, options);
    
    locate=y_opt;
end

xcor = linspace(5, 500, 100);
ycor = linspace(5, 500, 100);
z = 35;
array3D = zeros(length(xcor), length(ycor), 1);

for i = 1:length(xcor)
    for j = 1:length(ycor)
        s = [xcor(i), ycor(j), z];
        
        locate = find(s);
        % disp(s);
        % disp(S_est');
        array3D(i, j, :) = abs(norm(locate' - s));
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
zlabel('Error in location');
title('Localization percentage Error in 3D Space');
colorbar;



%% localization with the help of angle of arrival 

% assume actual source locations are given as this, we will estimate range of source with help of azimuthal, elevation angles and TDOA of sound signal for 2 microphones  
range = 60; % distance from origin
theta0 = 40; % polar angle in degrees
phi0 = 60; % azimuthal angle in degrees

% Convert angles from degrees to radians
theta0 = deg2rad(theta0);
phi0 = deg2rad(phi0);
mic1=[0, 0, 0];
mic2=[-1,0,0];

source= [range*sin(theta0) * cos(phi0),range*sin(theta0) * sin(phi0),range*cos(theta0)];
c = 343;
ts= 1/(64000*8);
Delta_t= round(((norm(source-mic1)-norm(source-mic2))/c)/ts)*ts;
direction_vector = [sin(theta0)*cos(phi0), sin(theta0)*sin(phi0), cos(theta0)];

% Distance difference calculation
distance_difference = Delta_t * c;

% Calculate r using the distance difference and geometry
% The sound travels different distances to each microphone

% Solve for r such that |r * direction_vector - mic1| - |r * direction_vector - mic2| = distance_difference
syms r;
eq = norm(r * direction_vector - mic1) - norm(r * direction_vector - mic2) == distance_difference;
sol = vpasolve(eq, r);

% Display the result
disp('The distance r is:');
disp(double(sol));

%%
% 
% Function to Estimate Distance (findrange):
% 
% The function findrange(delta_t, mic1, mic2, theta, phi) calculates the estimated distance (dist) to a source based on the time difference delta_t between two microphones (mic1 and mic2), and the angles theta (elevation) and phi (azimuth).
% It uses symbolic math to solve the equation that models the difference in distances from the source to each microphone, incorporating the direction vector defined by theta and phi.

% Simulation Parameters:
% 
% The code simulates a range of distances (range) for the sound source, varying from 5 to 500 meters.
% Microphones are positioned at [0, 0, 0] and [-1, 0, 0].
% A sampling period ts is derived from the sampling frequency to quantify the TDOA.
% The angles theta and phi are subject to random errors of up to ±2 degrees, simulating measurement inaccuracies.
% Range and Location Estimation:
% 
% For each distance in the range, the code calculates the source's position using both the actual and perturbed angles.
% The function findrange is called twice for each distance: once using the accurate angles (theta, phi) and once using the perturbed angles (theta_n, phi_n).
% The estimated range error is computed as the absolute difference between the estimated distance and the actual distance.
% The estimated location error is calculated as the Euclidean distance between the estimated and actual source positions.

function dist= findrange(delta_t, mic1, mic2,theta, phi)
    c=343;
    distance_difference = delta_t * c;
    
    direction_vector = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
    syms r;
    eq = norm(r * direction_vector - mic1) - norm(r * direction_vector - mic2) == distance_difference;
    sol = vpasolve(eq, r);
    dist= sol;
    
end


range= linspace(5,500,100);

mic1=[0, 0, 0];
mic2=[-1,0,0];
c = 343;
ts= 1/(64000*8*4);

estimated_range = zeros(length(range),1);
estimated_loaction_error= zeros(length(range),1);

for i=1:length(range)

    theta  = 50 ;
    phi= 60;

    theta_n =theta +(-2 + (4).*rand(1,1)); %there is random error in estimated theta of 2 degree
    phi_n= phi+(-2 + (4).*rand(1,1));  %there is random error in estimated phi of 2 degree
    theta_n= deg2rad(theta_n);
    phi_n= deg2rad(phi_n);
    theta = deg2rad(theta);
    phi = deg2rad(phi);
    % estimating the range when there is no error in theta and phi
    source= [range(i)*sin(theta) * cos(phi),range(i)*sin(theta) * sin(phi),range(i)*cos(theta)];
    Delta_t= round(((norm(source-mic1)-norm(source-mic2))/c)/ts)*ts;
    dist1= findrange(Delta_t,mic1,mic2,theta,phi);
    source_est= [dist1*sin(theta) * cos(phi),dist1*sin(theta) * sin(phi),dist1*cos(theta)];

    % estimating the range when there is error in theta and phi
    source_n=[range(i)*sin(theta_n) * cos(phi_n),range(i)*sin(theta_n) * sin(phi_n),range(i)*cos(theta_n)];
    Delta_t_n= round(((norm(source_n-mic1)-norm(source_n-mic2))/c)/ts)*ts;
    dist2=findrange(Delta_t_n,mic1,mic2,theta_n,phi_n);
    source_est_n= [dist2*sin(theta_n) * cos(phi_n),dist2*sin(theta_n) * sin(phi_n),dist2*cos(theta_n)];
    
    estimated_range(i,:)= abs(dist1-range(i));
    
    estimated_loaction_error(i,:)=norm(source_est_n-source);

end

figure 
stem(range,estimated_range);
xlabel(['actual range for theta = ',num2str(rad2deg(theta)),' and phi = ',num2str(rad2deg(phi))]);
ylabel('absolute estimated range error');

figure
stem(range,estimated_loaction_error);
xlabel(['actual range for theta with error upto 2 degree = ',num2str(rad2deg(theta)),' and phi = ',num2str(rad2deg(phi))])
ylabel('estimated error in location');




%%
%% comparing 6 microphone system of finding localization with AVS method of localization


function dist= findrange1(delta_t, mic1, mic2,theta, phi)
    c=343;
    distance_difference = delta_t * c;

    direction_vector = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
    syms r;
    eq = norm(r * direction_vector - mic1) - norm(r * direction_vector - mic2) == distance_difference;
    sol = vpasolve(eq, r);
    dist= sol;

end

function [r, theta, phi] = cartesianToSpherical(x, y, z)
    % Calculate radius
    r = sqrt(x.^2 + y.^2 + z.^2);

    % Calculate polar angle theta (angle from the z-axis)
    theta = acos(z ./ r);

    % Calculate azimuthal angle phi (angle from the x-axis in the xy-plane)
    phi = atan2(y, x); % atan2 handles the signs of x and y correctly

    % Display the results
    % fprintf('r = %.4f\n', r);
    % fprintf('theta = %.4f radians\n', theta);
    % fprintf('phi = %.4f radians\n', phi);
end

%%

mic1=[0, 0, 0];
mic2=[-1,0,0];
c = 343;
ts= 1/(64000*8*4);


s=[10, 77, 35];
[r, theta, phi]= cartesianToSpherical(s(1),s(2),s(3));

Delta_t= round(((norm(s-mic1)-norm(s-mic2))/c)/ts)*ts;


dist1= findrange1(Delta_t,mic1,mic2,theta,phi);


source_est = [dist1*sin(theta) * cos(phi),dist1*sin(theta) * sin(phi),dist1*cos(theta)];



error= abs(norm(source_est-s));

fprintf('error = %.4f\n', error)



%%



xcor = linspace(10, 100, 19);
ycor = linspace(10, 100, 19);
z = 35;
array3D = zeros(length(xcor), length(ycor), 1);

for i = 1:length(xcor)
    for j = 1:length(ycor)
        s = [xcor(i), ycor(j), z];
        [r, theta, phi]= cartesianToSpherical(s(1),s(2),s(3)); % angles in radian

        Delta_t= round(((norm(s-mic1)-norm(s-mic2))/c)/ts)*ts;

        dist1= findrange1(Delta_t,mic1,mic2,theta,phi);

        source_est = [dist1*sin(theta) * cos(phi),dist1*sin(theta) * sin(phi),dist1*cos(theta)];

        error1= abs(norm(source_est-s));

        array3D(i, j, :) =error1 ;
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
zlabel('Error in location');
title('Localization percentage Error in 3D Space');
colorbar;