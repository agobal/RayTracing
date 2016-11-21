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

x_laser = 0.0001;
y_laser = 0.0001;
z_laser = 0.0003;

intensity = 1;
alpha = 0.4;

while intensity > 0.01
    phi = pi;
    theta = pi/4;
    intersect = [];
    for i = 1:num_particles
        a = cos(theta)*sin(phi)/cos(phi);
        b = sin(theta)*sin(phi)/cos(phi);
        x0 = x(1, i);
        y0 = y(1, i);
        z0 = z(1, i);
        r0 = r(1, i);
        aa = x_laser - a*z_laser - x0;
        bb = y_laser - b*z_laser - y0;
        A = a^2 + b^2 + 1;
        B = (2*a*aa + 2*b*bb - 2*z0);
        C = aa^2 + bb^2 + z0^2 - r0^2;
        delta = B^2 - 4*A*C;
        if delta <= 0
            continue;
        else
            z1 = (-B + sqrt(delta))/(2*A);
            z2 = (-B - sqrt(delta))/(2*A);
            x1 = A*z1;
            y1 = B*z1;
            x2 = A*z2;
            y2 = B*z2;
            distance1 = sqrt((x1 - x_laser)^2 + (y1 - y_laser)^2 + (z1 - z_laser)^2);
            distance2 = sqrt((x2 - x_laser)^2 + (y2 - y_laser)^2 + (z2 - z_laser)^2);
            if distance1 < distance2
                intersect = [intersect; x1 y1 z1 i distance1];
            else
                intersect = [intersect; x2 y2 z2 i distance2];
            end
        end
    end
    [val indx] = min(intersect(:, 5));
    x_laser = intersect(indx, 1);
    y_laser = intersect(indx, 2);
    z_laser = intersect(indx, 3);
    thetaphiok = 0;
    while thetaphiok == 0
        theta = rand*2*pi;
        phi = rand*phi;
        coss = (x_laser - x(1, indx))*cos(theta)*sin(phi) + (y_laser - y(1, indx))*sin(theta)*sin(phi) + (z_laser - z(1, indx))*cos(phi);
        if coss <= 0
            thetaphiok = 1;
        else
            thetaphiok = 0;
        end
    end
    intensity = intensity*0.4
end