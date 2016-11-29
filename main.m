% %% 2. Create the two-scale particle packing
% 
% % Asking the user for inputs on bounding box
% 
% prompt = {'Model length (mm):', 'Model width (mm):', 'Model height (mm):', 'Layer thickness (mm):'};
% dlg_title = 'Model dimensions';
% num_lines = 1;
% defaultans = {'0.15', '0.15', '0.15', '0.15'};
% answer = inputdlg(prompt, dlg_title, num_lines, defaultans);
% 
% % Dimensions of the packing box
% x_container = str2double(char(answer(1)))/1000;
% y_container = str2double(char(answer(2)))/1000;
% z_container = str2double(char(answer(3)))/1000;
% 
% % Layer thickness
% z_layer = str2double(char(answer(4)))/1000;
% 
% % 2.1. Creating the smaller scale (particle level) packing 
% 
% % Asking the user for inputs on powder particle properties
% 
% prompt = {'Average particle radius (mm):', 'Standard deviation:', 'Packing fraction:'};
% dlg_title = 'Powder properties';
% num_lines = 1;
% defaultans = {'0.0075', '0.003', '0.64'};
% answer = inputdlg(prompt, dlg_title, num_lines, defaultans);
% 
% average_particle_radius = str2double(char(answer(1)))/1000;
% stddev = str2double(char(answer(2)))/1000;
% packing_density = str2double(char(answer(3)));
% 
% 
% 
% disp('Creating the particle packing');
% [x, y, z, r, num_particles] = PackingGenerator(x_container, y_container, z_container, average_particle_radius, stddev, packing_density);


%% Shooting rays inside the packing from random locations on the top layer

x1 = 0.0001;
y1 = 0.0001;
z1 = 0.0001;

x2 = 0.0001;
y2 = 0.0001;
z2 = 0.0000;

intensity = 1;
alpha = 0.4;

for j = 1:1000000
    j
    intensity = 1;
    while intensity > 0.01
        intersect = [];
        for i = 1:num_particles
            x0 = x(1, i);
            y0 = y(1, i);
            z0 = z(1, i);
            r0 = r(1, i);
            % Find a, b and c for solving the equation
            a = (x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2;
            b = 2*((x2 - x1)*(x1 - x0) + (y2 - y1)*(y1 - y0) + (z2 - z1)*(z1 - z0));
            c = x0^2 + y0^2 + z0^2 + x1^2 + y1^2 + z1^2 - 2*(x0*x1 + y0*y1 + z0*z1) - r0^2;
            delta = b^2 - 4*a*c;
            if delta <= 0
                continue;
            else
                t1 = (-b + sqrt(delta))/(2*a);
                t2 = (-b - sqrt(delta))/(2*a);
                x1p = x1 + (x2 - x1)*t1;
                y1p = y1 + (y2 - y1)*t1;
                z1p = z1 + (z2 - z1)*t1;
                x2p = x1 + (x2 - x1)*t2;
                y2p = y1 + (y2 - y1)*t2;
                z2p = z1 + (z2 - z1)*t2;
                distance1 = sqrt((x1p - x1)^2 + (y1p - y1)^2 + (z1p - z1)^2);
                distance2 = sqrt((x2p - x1)^2 + (y2p - y1)^2 + (z2p - z1)^2);
                if distance1 < distance2
                    intersect = [intersect; x1p y1p z1p i distance1];
                else
                    intersect = [intersect; x2p y2p z2p i distance2];
                end
            end
        end
        [val indx] = min(intersect(:, 5));
        x1 = intersect(indx, 1);
        y1 = intersect(indx, 2);
        z1 = intersect(indx, 3);
        theta = rand*2*pi;
        phi = rand*pi;
        x2 = x0 + r0*sin(phi)*cos(theta);
        y2 = y0 + r0*sin(phi)*sin(theta);
        z2 = z0 + r0*cos(phi);
        intensity = intensity*0.4;
    end
end