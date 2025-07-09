close all
clearvars


%Q3/4/5
kb = 1;

lambda = 2 * pi / kb;   %absolute locations
h = lambda / 20;

%expressing locations in terms of h:
Nx = 21;       %[lambda, lambda]
Ny = 21;
source_location = reshape([11, 201], [1, 1, 2]);    %[lambda/2, 10*lambda]
x = 1:Nx;
y = 1:Ny;   % == Nx in this case!

[X,Y] = meshgrid(x,y);
N_tot = numel(X); %number of grid elements



% Q6
%need some sort of distance metric:
grid_tensor = reshape([X, Y], [Nx, Ny, 2] ); 
displacement_tensor_to_source = source_location - grid_tensor;     
distances_to_source = h * vecnorm(displacement_tensor_to_source, 2, 3);    %making use of tensor stuff :)

%hankel computation
u_inc = -1j / 4 * hankel(kb * distances_to_source);





run("genimage.m")



lambda_tikh = 1e-9;


% Q10/11
M_per_array = 25; %per array! more in total


rec_locations_x_original = linspace(-19, 41, M_per_array); %origin is in (1,1), so -lambda == -19 in coords
rec_locations_y_original = linspace(31, 31, M_per_array);


%add more receivers:
rec_locations_x_west = linspace(-19, -19, M_per_array);
rec_locations_y_west = linspace(31, -11, M_per_array);

rec_locations_x_east = linspace(40, 40, M_per_array);
rec_locations_y_east = linspace(31, -11, M_per_array);

rec_locations_x_north = linspace(-19, 41, M_per_array);
rec_locations_y_north = linspace(-11, -11, M_per_array);

rec_locations_x = [rec_locations_x_original rec_locations_x_west rec_locations_x_east rec_locations_x_north];
rec_locations_y = [rec_locations_y_original rec_locations_y_west rec_locations_y_east rec_locations_y_north];

[~, M] = size(rec_locations_x);
rec_locations = zeros(1, 1, 2, M);

rec_locations(:,:,1,:) = rec_locations_x; %to prep for tensor stuff
rec_locations(:,:,2,:) = rec_locations_y;



x_contrast = image(:);




displacement_tensor_to_rec = grid_tensor - rec_locations; %(Nx, Nx , 2, M); %[grid_matrix x 2D x num_receivers]
distances_to_rec = squeeze(h * vecnorm(displacement_tensor_to_rec, 2, 3));   

distances_to_source_flattened = distances_to_source(:);
distances_to_rec_flattened = reshape(distances_to_rec, [N_tot, M]);




%vertically stack them!
kb_min = 1;
kb_max = 10;
kb_points = 20;



A = [];
for kb = linspace(kb_min, kb_max, kb_points)
    A_kb = zeros(M, N_tot);
    for i = 1:M
        for j = 1:N_tot
            A_kb(i,j) = -1/16 * h^2 * hankel(kb *  distances_to_rec_flattened(j, i)) * hankel(kb * distances_to_source_flattened(j));
    
        end
    end

    A = vertcat(A, A_kb);
end




s = svd(A);

% Plot them on a log scale
figure;
semilogy(s, 'bo-', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Singular Value Index', 'FontSize', 14);
ylabel('Magnitude (log scale)', 'FontSize', 14);
title('Singular Values of Matrix A', 'FontSize', 16);
grid on;
set(gca, 'FontSize', 20);

    %note: what happens to the singular values here??? for in report





sigma = 0.00001;   %0.0001; % or 0.0001
u_sc_noiseless = A*x_contrast;  %this is what we measure
noise = sigma*(randn(size(u_sc_noiseless))  +1j*randn(size(u_sc_noiseless))   ) / sqrt(2);
u_sc = u_sc_noiseless + noise;




% pinv based reconstructions
x_estimate_pinv = pinv(A) * u_sc;
image_pinv_rec = reshape(x_estimate_pinv, [Nx Ny]);

AH = A';
x_estimate_tikh = (AH * A + lambda_tikh * eye(size(A, 2))) \ (AH * u_sc);
image_tikh_rec = reshape(x_estimate_tikh, [Nx Ny]);




figure
imagesc(abs(image_pinv_rec))
colorbar()
title('Pinv Reconstruction')
colormap(gray);
axis equal tight
set(gca, 'FontSize', 20);  % axes tick labels

figure
imagesc(abs(image_tikh_rec))
colorbar()
title(['Tikhonov Regularized Reconstruction, \lambda = ' num2str(lambda_tikh)])
colormap(gray);
axis equal tight
set(gca, 'FontSize', 20);  % axes tick labels






%adding rotations:







