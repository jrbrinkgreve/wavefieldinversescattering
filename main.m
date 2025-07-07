%main file for 2D inverse scattering problem
close all
clearvars
%Q1:
%the hankel function is a linear combination of bessel functions. As bessel
%functions are the solutions to the wave equation in cylindrical
%coordinates, it does not surprise me this is the case.





%Q2
%see notebook!



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

%plotting:
figure;
% Plot real part
subplot(1,2,1);
imagesc(real(u_inc));
colorbar;
title('Real Part of u\_inc');
axis equal tight;

% Plot imaginary part
subplot(1,2,2);
imagesc(imag(u_inc));
colorbar;
title('Imaginary Part of u\_inc');
axis equal tight;




% Q7
%changing parameters and observations:

%case source closer: the field initially becomes stronger and the phase
%shifts a bit. When the source is really close, you can see the wavefront
%becoming spherical! (include plots in report!)

%changing kb to 2: 
%as lambda is defined as 2pi/kb, and the displacements are defined in terms
%of lambda there is nothing happening when kb is changed. In the hankel
%expression it can be interpreted as if the source would be further away,
%but the sources move physically closer as the wavelength shrinks. So this
%effect is cancelled. Would the physical distances remain the same, then it
%would be like the field is from further away with larger kb.








% Q8 / 9
%get image
run("genimage.m")



tikh_images = zeros(Nx,Ny, 4);
pinv_images = zeros(Nx,Ny, 4);

for k = 0:0     %deg = [0 90 180 270] rotations
image = rot90(image_original, k);

% Q10/11
M_per_array = 30; %per array! more in total


rec_locations_x_original = linspace(-19, 41, M_per_array); %origin is in (1,1), so -lambda == -19 in coords
rec_locations_y_original = linspace(31, 31, M_per_array);

%add more receivers:
rec_locations_x_west = linspace(-19, -19, M_per_array);
rec_locations_y_west = linspace(20, -20, M_per_array);

rec_locations_x_east = linspace(40, 40, M_per_array);
rec_locations_y_east = linspace(20, -20, M_per_array);

rec_locations_x_north = linspace(-19, 41, M_per_array);
rec_locations_y_north = linspace(-31, -31, M_per_array);


rec_locations = zeros(1, 1, 2, 4*M_per_array);
rec_locations_x = [rec_locations_x_original rec_locations_x_west rec_locations_x_east rec_locations_x_north];
rec_locations_y = [rec_locations_y_original rec_locations_y_west rec_locations_y_east rec_locations_y_north];

[~, M] = size(rec_locations_x);

rec_locations(:,:,1,:) = rec_locations_x; %to prep for tensor stuff
rec_locations(:,:,2,:) = rec_locations_y;



% Q12
%see notebook:

% Q13
x_contrast = image(:); %no need for reshape to vectorise(:)


% Q14
%also need distance from image dom to  rec: for every point in dom, to
%every mic = M*N_tot combinations: and also 2-d, so need another dim for
%storing displacements and distances

displacement_tensor_to_rec = grid_tensor - rec_locations; %(Nx, Nx , 2, M); %[grid_matrix x 2D x num_receivers]
distances_to_rec = squeeze(h * vecnorm(displacement_tensor_to_rec, 2, 3));   

distances_to_source_flattened = distances_to_source(:);
distances_to_rec_flattened = reshape(distances_to_rec, [N_tot, M]);

A = zeros(M, N_tot);
for i = 1:M
    for j = 1:N_tot
        A(i,j) = -1/16 * h^2 * hankel(kb *  distances_to_rec_flattened(j, i)) * hankel(distances_to_source_flattened(j));

    end
end
%writing convolution integral as a convolution sum, and then writing the
%conv sum as a matrix multiplication!


% Q15
singular_values = svd(A, 'econ');
singular_values(1:min(M, 20));
    %note: what happens to the singular values here??? for in report


% Q16
u_sc = A*x_contrast;  %this is what we measure





% Q17 / 18
%invert the system: 



x_estimate_pinv = pinv(A) * u_sc;
image_pinv_rec = reshape(x_estimate_pinv, [Nx Ny]);
pinv_images(:,:,k+1) = rot90(image_pinv_rec, -k);

%{
figure

imagesc((abs(image_pinv_rec)))
colorbar()
title(['Pinv Reconstruction'] )
axis equal tight
%}


% Q19
%more receivers: add plots in report!

%note: know there might be noise on top due to the code required for the
%next section



%Q20: add noise
%do analysis... seems pretty stable!





% Q21: try and improve everything!
%notes: what to try:
% - add receivers around the object instead of just under
% - add multiple source locations / rotate the measurement orientation and
%           add that to the measurement to illuminate from multiple
%           directions!
% - use another reconstruction method instead of simple rec...

%i think this question is what determines the grade of the assignment!



%current state of the code: Tikhonov regulatization done + more sensors
%around the object





lambda = 1e-7; % regularization parameter
AH = A';
x_estimate_tikh = (AH * A + lambda * eye(size(A, 2))) \ (AH * u_sc);
image_tikh_rec = reshape(x_estimate_tikh, [Nx Ny]);

tikh_images(:,:, k+1) = rot90(image_tikh_rec, -k);
avg_estimate = mean(tikh_images, 3);

end





%note: first check neumann convergence criteria
%note: could define object operator, and do neumann iterations




%for more advanced schemes we need the object operator Gd:
%-------------------



%just a bunch of tensor abuse below
inv_dom_locations_x = grid_tensor(:,:,1);
inv_dom_locations_x = inv_dom_locations_x(:);

inv_dom_locations_y = grid_tensor(:,:,2);
inv_dom_locations_y = inv_dom_locations_y(:);

inv_dom_locations = zeros(1, 1, 2, N_tot);
inv_dom_locations(:,:,1,:) = inv_dom_locations_x;
inv_dom_locations(:,:,2,:) = inv_dom_locations_y;

displacement_tensor_inv_to_inv = grid_tensor - inv_dom_locations;
distances_inside_inv = squeeze(h * vecnorm(displacement_tensor_inv_to_inv, 2, 3)); 
distances_inside_inv_flattened = reshape(distances_inside_inv, [N_tot, N_tot]);



%for more advanced schemes we need the object operator:
B = zeros(N_tot, N_tot);
for i = 1:N_tot
    for j = 1:N_tot
        B(i,j) = -1j/4 * h^2 * hankel(kb *  distances_inside_inv_flattened(j, i));

    end
end

for i = 1:N_tot
    B(i,i) = 0; %to get rid of hankel(0) == NaN
end







%init
max_iter = 500;
chi = zeros(Nx);
chi = chi(:);
u_tot = u_inc(:); 


%iterate:
for k = 1:max_iter



%updated A matrix as u changes
A = (-1/16) * h^2 * (u_tot(:) .* hankel(kb * distances_to_rec_flattened)).';

%(prev unvectorized A calc)

%{
A = zeros(M, N_tot);
for i = 1:M
    for j = 1:N_tot
        A(i,j) = -1/16 * h^2 * hankel(kb *  distances_to_rec_flattened(j, i)) * u_tot(j);

    end
end
%}


 % predicted scattered field
u_sc_est = A * chi;          
residue = u_sc - u_sc_est;        % data residual, u_sc is 'true' measurement

%solve system
AH = A';
delta_chi = (AH*A + lambda*eye(N_tot)) \ (AH * residue);

%update contrast
chi = chi + delta_chi;

%new estimate of field
scatter_source = chi .* u_tot;        % contrast × field
u_tot = u_inc(:) + B * scatter_source;



%convergence tracking
res_norm = norm(residue) / norm(u_sc);
if mod(k,100) == 0
fprintf('Iteration %d: relative residual = %.4e\n', k, res_norm);
end
if res_norm < 1e-9
    break;
end
end

%form image
reconstructed_image = reshape(chi, [Nx, Ny]);

%display it
figure;
imagesc(abs(reconstructed_image));
title('DWBI Reconstructed Image');
colorbar;
axis equal tight;
















figure
imagesc(abs(pinv_images(:,:,1)))
colorbar()
title('Pinv Reconstruction')
axis equal tight

figure
imagesc(abs(tikh_images(:,:,1)))
colorbar()
title(['Tikhonov Regularized Reconstruction, \lambda = ' num2str(lambda)])
axis equal tight


figure
imagesc(abs(avg_estimate))
colorbar()
title(['Tikhonov + rotated illumination, \lambda = ' num2str(lambda)])
axis equal tight


