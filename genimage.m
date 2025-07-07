%genimage file
image = zeros(Nx, Ny);









%ellipse, or a potato :)
%{

%center of ellipse
cx = 11; cy = 11;

% Ellipse radii
a = 8; % major axis
b = 4; % minor axis

% Rotation angle in radians (e.g. 30 degrees)
theta = pi/6; 

% Loop over every pixel
for x = 1:Nx
    for y = 1:Ny
        % Shift coordinates to center
        xp = x - cx;
        yp = y - cy;
        
        % Rotate coordinates by -theta
        Xp =  cos(theta)*xp + sin(theta)*yp;
        Yp = -sin(theta)*xp + cos(theta)*yp;
        
        % Check ellipse equation
        val = (Xp^2)/(a^2) + (Yp^2)/(b^2);
        
        if val <= 1
            image(y,x) = 1;
        end
    end
end
%}











%ball:
%{

% Define center
cx = 11; cy = 11; % center of the matrix

% Fill circle with radius 5
for x = 1:21
    for y = 1:21
        r = sqrt((x - cx)^2 + (y - cy)^2);
        if r <= 5
            image(y, x) = 1;
        end
    end
end

image = image * 0.01; % to reduce contrast
%}









%square:
%{
side = 9;

% Define start and end indices
start_idx = floor((21 - side)/2) + 1;
end_idx = start_idx + side - 1;

% Fill the square
image(start_idx:end_idx, start_idx:end_idx) = 1;

%}






%half-circle:
% Parameters
radius = 6;
center = ceil(Nx / 2);  % Center at (11,11) if Nx = 21
tolerance = 0.5;  % Thickness of the arc



% Draw hollow upper half-circle
for x = 1:Nx
    for y = 1:Nx
        dx = x - center;
        dy = y - center;
        distance = sqrt(dx^2 + dy^2);

        % Check for hollow arc in upper half
        if abs(distance - radius) <= tolerance && dy <= 0
            image(y, x) = 1;
        end
    end
end






%--------------------------


%in case we also want noisy images, but shouldnt matter too much

image_no_noise = image;

image_no_noise = image_no_noise * 0.1; %for Born
sigma = 0;%0.01;
noise = sigma  * randn(size(image));

image_original = image_no_noise + noise;

%Q9 display


figure;
% Plot real part
imagesc(image_no_noise);
colormap(gray);
colorbar
axis equal tight;
title('Object')

