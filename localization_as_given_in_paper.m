global c s1 s2 s3 s4 s5 s6 s7 s8 Ts;
c = 343;
%%
% Define the microphone positions for both acoustic arrays
s1 = [0, 0, 0];
s2 = [0, 1, 0];
s3 = [1, 0, 0];
s4 = [0, 0, 1];
s5 = [0, 0, 5];
s6 = [0, 0, 6];
s7 = [-1, 0, 5];
s8 = [0, -1, 5];
s =[40,30,110]; % original source location

Ts=1/(64000*8);

% Microphone positions array
mic_positions = [s1; s2; s3; s4; s5; s6; s7; s8];

% TDOA= [(norm(s-s2)-norm(s-s1))/c,(norm(s-s3)-norm(s-s1))/c,(norm(s-s4)-norm(s-s1))/c,(norm(s-s6)-norm(s-s5))/c,(norm(s-s7)-norm(s-s5))/c,(norm(s-s8)-norm(s-s5))/c ];
TDOA= [round(((norm(s-s2)-norm(s-s1))/c)/Ts)*Ts,round(((norm(s-s3)-norm(s-s1))/c)/Ts)*Ts,round(((norm(s-s4)-norm(s-s1))/c)/Ts)*Ts,round(((norm(s-s6)-norm(s-s5))/c)/Ts)*Ts,round(((norm(s-s7)-norm(s-s5))/c)/Ts)*Ts,round(((norm(s-s8)-norm(s-s5))/c)/Ts)*Ts];
disp(TDOA);

% Example TDOA values (replace these with actual TDOA calculations)
% TDOA = [t21, t31, t41, t65, t75, t85]; % Replace with actual TDOA data


% Form the A matrix and B vector for least squares calculation
% A = 2 * [
%     s2 - s1, c*t21, 0;
%     s3 - s1, c*t31, 0;
%     s4 - s1, c*t41, 0;
%     s6 - s5, 0, c*t65;
%     s7 - s5, 0, c*t75;
%     s8 - s5, 0, c*t85
% ];
% 
% B = [
%     norm(s2)^2 - norm(s1)^2 - (c*t21)^2;
%     norm(s3)^2 - norm(s1)^2 - (c*t31)^2;
%     norm(s4)^2 - norm(s1)^2 - (c*t41)^2;
%     norm(s6)^2 - norm(s5)^2 - (c*t65)^2;
%     norm(s7)^2 - norm(s5)^2 - (c*t75)^2;
%     norm(s8)^2 - norm(s5)^2 - (c*t85)^2
% ];
% 
% % Solve for the drone's position using least squares
% XLos = (A' * A) \ (A' * B);
% 
% % Extract the position coordinates
% drone_position = XLos(1:3);
% fprintf('Estimated Drone Position: X = %.2f, Y = %.2f, Z = %.2f\n', drone_position);

% Optional: Refine the position iteratively (if necessary)
% Implement iterative refinement here if required based on initial XLos


%%
% Perform the localization calculations
A = 2 * [
    s2 - s1, c*TDOA(1), 0;
    s3 - s1, c*TDOA(2), 0;
    s4 - s1, c*TDOA(3), 0;
    s6 - s5, 0, c*TDOA(4);
    s7 - s5, 0, c*TDOA(5);
    s8 - s5, 0, c*TDOA(6)
];

B = [
    norm(s2)^2 - norm(s1)^2 - (c*TDOA(1))^2;
    norm(s3)^2 - norm(s1)^2 - (c*TDOA(2))^2;
    norm(s4)^2 - norm(s1)^2 - (c*TDOA(3))^2;
    norm(s6)^2 - norm(s5)^2 - (c*TDOA(4))^2;
    norm(s7)^2 - norm(s5)^2 - (c*TDOA(5))^2;
    norm(s8)^2 - norm(s5)^2 - (c*TDOA(6))^2
];

% Solve for the drone's position using least squares
XLos = (A' * A) \ (A' * B);

% Extract the position coordinates
drone_position = XLos(1:3);
fprintf('Estimated Drone Position: X = %.2f, Y = %.2f, Z = %.2f\n', drone_position);


%%
function TDOAarray=findTDOA(s)
    global c s1 s2 s3 s4 s5 s6 s7 s8 Ts;
    TDOAarray = [
        round(((norm(s - s2) - norm(s - s1)) / c) / Ts) * Ts, ...
        round(((norm(s - s3) - norm(s - s1)) / c) / Ts) * Ts, ...
        round(((norm(s - s4) - norm(s - s1)) / c) / Ts) * Ts, ...
        round(((norm(s - s6) - norm(s - s5)) / c) / Ts) * Ts, ...
        round(((norm(s - s7) - norm(s - s5)) / c) / Ts) * Ts, ...
        round(((norm(s - s8) - norm(s - s5)) / c) / Ts) * Ts
    ];
end


function location= findlocation(tdoa)
    global c s1 s2 s3 s4 s5 s6 s7 s8 ;
    
    A = 2 * [
        s2 - s1, c*tdoa(1), 0;
        s3 - s1, c*tdoa(2), 0;
        s4 - s1, c*tdoa(3), 0;
        s6 - s5, 0, c*tdoa(4);
        s7 - s5, 0, c*tdoa(5);
        s8 - s5, 0, c*tdoa(6)
    ];
    
    B = [
        norm(s2)^2 - norm(s1)^2 - (c*tdoa(1))^2;
        norm(s3)^2 - norm(s1)^2 - (c*tdoa(2))^2;
        norm(s4)^2 - norm(s1)^2 - (c*tdoa(3))^2;
        norm(s6)^2 - norm(s5)^2 - (c*tdoa(4))^2;
        norm(s7)^2 - norm(s5)^2 - (c*tdoa(5))^2;
        norm(s8)^2 - norm(s5)^2 - (c*tdoa(6))^2
    ];
    
    % Solve for the drone's position using least squares
    XLos = (A' * A) \ (A' * B);
    
    % Extract the position coordinates
    location = XLos(1:3);
    
end

x=linspace(10,110,11);
y=linspace(10,110,11);
z=linspace(10,110,11);

array3D = zeros(length(x), length(y), length(z),3);

for i=1:length(x)
    for j=1:length(y)
        for k= 1:length(z)
            s= [x(i) ,y(j) ,z(k)];
            tdoa= findTDOA(s);
            array3D(i,j,k,:)=findlocation(tdoa);
        end
    end

end
%%



% Flatten the arrays for plotting
[X, Y, Z] = ndgrid(x, y, z);
X = X(:);
Y = Y(:);
Z = Z(:);
calcX = array3D(:, :, :, 1);
calcY = array3D(:, :, :, 2);
calcZ = array3D(:, :, :, 3);
calcX = calcX(:);
calcY = calcY(:);
calcZ = calcZ(:);

% Create a 3D scatter plot of the original points
figure;
scatter3(X, Y, Z, 'b');
hold on;

% Create a 3D scatter plot of the calculated locations
scatter3(calcX, calcY, calcZ, 'r');

% Set labels and legend
xlabel('X');
ylabel('Y');
zlabel('Z');
legend('Original Points', 'Calculated Locations');
title('3D Visualization of Original Points and Calculated Locations');
grid on;
hold off;


%%

%when reference for 1st array is s1 and s5 for 2nd
% Define the microphone positions for both acoustic arrays
s1 = [0, 0, 1];
s2 = [0, 0.3, 0];
s3 = [0.3, 0, 0];
s4 = [0, 0, 0.3];
s5 = [0, 0, 0];
s6 = [0, 0, 1.3];
s7 = [-0.3, 0, 1];
s8 = [0, -0.3, 1];
s =[60,70,100];

Ts=1/(64000*8);

% Microphone positions array
mic_positions = [s1; s2; s3; s4; s5; s6; s7; s8];

%TDOA= [(norm(s-s2)-norm(s-s5))/c,(norm(s-s3)-norm(s-s5))/c,(norm(s-s4)-norm(s-s5))/c,(norm(s-s6)-norm(s-s1))/c,(norm(s-s7)-norm(s-s1))/c,(norm(s-s8)-norm(s-s1))/c ];
TDOA= [round(((norm(s-s2)-norm(s-s1))/c)/Ts)*Ts,round(((norm(s-s3)-norm(s-s1))/c)/Ts)*Ts,round(((norm(s-s4)-norm(s-s1))/c)/Ts)*Ts,round(((norm(s-s6)-norm(s-s5))/c)/Ts)*Ts,round(((norm(s-s7)-norm(s-s5))/c)/Ts)*Ts,round(((norm(s-s8)-norm(s-s5))/c)/Ts)*Ts];
disp(TDOA);

% Perform the localization calculations
A = 2 * [
    s2 - s1, c*TDOA(1), 0;
    s3 - s1, c*TDOA(2), 0;
    s4 - s1, c*TDOA(3), 0;
    s6 - s5, 0, c*TDOA(4);
    s7 - s5, 0, c*TDOA(5);
    s8 - s5, 0, c*TDOA(6)
];

B = [
    norm(s2)^2 - norm(s1)^2 - (c*TDOA(1))^2;
    norm(s3)^2 - norm(s1)^2 - (c*TDOA(2))^2;
    norm(s4)^2 - norm(s1)^2 - (c*TDOA(3))^2;
    norm(s6)^2 - norm(s5)^2 - (c*TDOA(4))^2;
    norm(s7)^2 - norm(s5)^2 - (c*TDOA(5))^2;
    norm(s8)^2 - norm(s5)^2 - (c*TDOA(6))^2
];

% Solve for the drone's position using least squares
XLos = (A' * A) \ (A' * B);

% Extract the position coordinates
drone_position = XLos(1:3);
fprintf('Estimated Drone Position: X = %.2f, Y = %.2f, Z = %.2f\n', drone_position);


%%
%whens s5 is reference for first array and s1 is reference for 2nd array


% Define the microphone positions for both acoustic arrays
s1 = [0, 0, 0];
s2 = [0, 0.3, 0];
s3 = [0.3, 0, 0];
s4 = [0, 0, 0.3];
s5 = [0, 0, 2];
s6 = [0, 0, 2.3];
s7 = [-0.3, 0, 2];
s8 = [0, -0.3, 2];
s =[40,60,20];

Ts=1/(64000*8);

% Microphone positions array
mic_positions = [s1; s2; s3; s4; s5; s6; s7; s8];

%TDOA= [(norm(s-s2)-norm(s-s5))/c,(norm(s-s3)-norm(s-s5))/c,(norm(s-s4)-norm(s-s5))/c,(norm(s-s6)-norm(s-s1))/c,(norm(s-s7)-norm(s-s1))/c,(norm(s-s8)-norm(s-s1))/c ];
TDOA= [round(((norm(s-s2)-norm(s-s5))/c)/Ts)*Ts,round(((norm(s-s3)-norm(s-s5))/c)/Ts)*Ts,round(((norm(s-s4)-norm(s-s5))/c)/Ts)*Ts,round(((norm(s-s6)-norm(s-s1))/c)/Ts)*Ts,round(((norm(s-s7)-norm(s-s1))/c)/Ts)*Ts,round(((norm(s-s8)-norm(s-s1))/c)/Ts)*Ts];
disp(TDOA);

% Perform the localization calculations
A = 2 * [
    s2 - s5, c*TDOA(1), 0;
    s3 - s5, c*TDOA(2), 0;
    s4 - s5, c*TDOA(3), 0;
    s6 - s1, 0, c*TDOA(4);
    s7 - s1, 0, c*TDOA(5);
    s8 - s1, 0, c*TDOA(6)
];

B = [
    norm(s2)^2 - norm(s5)^2 - (c*TDOA(1))^2;
    norm(s3)^2 - norm(s5)^2 - (c*TDOA(2))^2;
    norm(s4)^2 - norm(s5)^2 - (c*TDOA(3))^2;
    norm(s6)^2 - norm(s1)^2 - (c*TDOA(4))^2;
    norm(s7)^2 - norm(s1)^2 - (c*TDOA(5))^2;
    norm(s8)^2 - norm(s1)^2 - (c*TDOA(6))^2
];

% Solve for the drone's position using least squares
XLos = (A' * A) \ (A' * B);

% Extract the position coordinates
drone_position = XLos(1:3);
fprintf('Estimated Drone Position: X = %.2f, Y = %.2f, Z = %.2f\n', drone_position);



%%

% 8 microphone system error analysis

%%


c = 343;

% Define the microphone positions for both acoustic arrays
s1 = [0, 0, 0];
s2 = [0, 0.5, 0];
s3 = [0.5, 0, 0];
s4 = [0, 0, 0.5];
s5 = [0, 0, 5];
s6 = [0, 0, 5.5 ];
s7 = [-0.5, 0, 5 ];
s8 = [0, -0.5, 5];
s =[40,50,20];

Ts=1/(64000*8);

% Microphone positions array
mic_positions = [s1; s2; s3; s4; s5; s6; s7; s8];

% TDOA= [(norm(s-s2)-norm(s-s1))/c,(norm(s-s3)-norm(s-s1))/c,(norm(s-s4)-norm(s-s1))/c,(norm(s-s6)-norm(s-s5))/c,(norm(s-s7)-norm(s-s5))/c,(norm(s-s8)-norm(s-s5))/c ];
TDOA= [round(((norm(s-s2)-norm(s-s1))/c)/Ts)*Ts,round(((norm(s-s3)-norm(s-s1))/c)/Ts)*Ts,round(((norm(s-s4)-norm(s-s1))/c)/Ts)*Ts,round(((norm(s-s6)-norm(s-s5))/c)/Ts)*Ts,round(((norm(s-s7)-norm(s-s5))/c)/Ts)*Ts,round(((norm(s-s8)-norm(s-s5))/c)/Ts)*Ts];
disp(TDOA);

%%
% Perform the localization calculations
A = 2 * [
    s2 - s1, c*TDOA(1), 0;
    s3 - s1, c*TDOA(2), 0;
    s4 - s1, c*TDOA(3), 0;
    s6 - s5, 0, c*TDOA(4);
    s7 - s5, 0, c*TDOA(5);
    s8 - s5, 0, c*TDOA(6)
];

B = [
    norm(s2)^2 - norm(s1)^2 - (c*TDOA(1))^2;
    norm(s3)^2 - norm(s1)^2 - (c*TDOA(2))^2;
    norm(s4)^2 - norm(s1)^2 - (c*TDOA(3))^2;
    norm(s6)^2 - norm(s5)^2 - (c*TDOA(4))^2;
    norm(s7)^2 - norm(s5)^2 - (c*TDOA(5))^2;
    norm(s8)^2 - norm(s5)^2 - (c*TDOA(6))^2
];

% Solve for the drone's position using least squares
XLos = (A' * A) \ (A' * B);

% Extract the position coordinates
drone_position = XLos(1:3);
fprintf('Estimated Drone Position: X = %.2f, Y = %.2f, Z = %.2f\n', drone_position);


%%
function TDOAarray=findTDOA2(source)
    global c s1 s2 s3 s4 s5 s6 s7 s8 Ts;
    TDOAarray = [
        round(((norm(source - s2) - norm(source - s1)) / c) / Ts) * Ts, ...
        round(((norm(source - s3) - norm(source - s1)) / c) / Ts) * Ts, ...
        round(((norm(source - s4) - norm(source - s1)) / c) / Ts) * Ts, ...
        round(((norm(source - s6) - norm(source - s5)) / c) / Ts) * Ts, ...
        round(((norm(source - s7) - norm(source - s5)) / c) / Ts) * Ts, ...
        round(((norm(source - s8) - norm(source - s5)) / c) / Ts) * Ts
    ];
end


function location= findlocation2(tdoa)
    global c s1 s2 s3 s4 s5 s6 s7 s8 ;
    
    A = 2 * [
        s2 - s1, c*tdoa(1), 0;
        s3 - s1, c*tdoa(2), 0;
        s4 - s1, c*tdoa(3), 0;
        s6 - s5, 0, c*tdoa(4);
        s7 - s5, 0, c*tdoa(5);
        s8 - s5, 0, c*tdoa(6)
    ];
    
    B = [
        norm(s2)^2 - norm(s1)^2 - (c*tdoa(1))^2;
        norm(s3)^2 - norm(s1)^2 - (c*tdoa(2))^2;
        norm(s4)^2 - norm(s1)^2 - (c*tdoa(3))^2;
        norm(s6)^2 - norm(s5)^2 - (c*tdoa(4))^2;
        norm(s7)^2 - norm(s5)^2 - (c*tdoa(5))^2;
        norm(s8)^2 - norm(s5)^2 - (c*tdoa(6))^2
    ];
    
    % Solve for the drone's position using least squares
    XLos = (A' * A) \ (A' * B);
    
    % Extract the position coordinates
    location = XLos(1:3);
    
end

x=linspace(5,500,100);
y=linspace(5,500,100);
z=40;

array3D = zeros(length(x), length(y),1);


for i=1:length(x)
    for j=1:length(y)
        s= [x(i) ,y(j) ,z];
        tdoa= findTDOA2(s);
        S_est=findlocation2(tdoa);
        % disp(s);
        % disp(S_est');
        array3D(i,j,:)=abs(norm(S_est'-s));
    end
end


%% Plot the data
[X, Y] = meshgrid(x, y);
Z = squeeze(array3D(:,:,1));

figure;
surf(X, Y, Z);
xlabel('X Position');
ylabel('Y Position');
zlabel('Error');
title('Localization Error in 3D Space');
colorbar;

