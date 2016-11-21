function [x, y, z, r, num_particles] = PackingGenerator(x_container, y_container, z_container, average_particle_radius, stddev, packing_density)
tic


% First step is to create the particle size distribution
V_particles = 0;
for i = 1:200000
    r(i, 1) = stddev*rand + average_particle_radius;
    V_particles = V_particles + (4/3)*pi*(r(i, 1)^3);
    if V_particles >= (x_container*y_container*z_container*packing_density)
        num_particles = i - 1;
        break;
    end
end
r = r(1:num_particles, 1)';


% Initialize particle positions
x = zeros(1, num_particles);
y = zeros(1, num_particles);
z = zeros(1, num_particles);
overlap = zeros(1, num_particles);

%Random positioning of the particles in the bounding box
for j = 1:num_particles
    x(1, j) = rand*(x_container - 2*r(1, j)) + r(1, j);
    y(1, j) = rand*(y_container - 2*r(1, j)) + r(1, j);
    z(1, j) = rand*(z_container - 2*r(1, j)) + r(1, j);
end


%Relocation of particles to reduce overlaps
for k = 1:1000 % Number of relocation steps
    k
    for i = 1:num_particles
        % Intermediate position of particles
        q = 0;
        xn = 0;
        yn = 0;
        zn = 0;
        for j = 1:num_particles
            if i ~= j
                overlap(1, j) = (r(1, i) + r(1, j) - sqrt((x(1, i) - x(1, j))^2 + (y(1, i) - y(1, j))^2 + (z(1, i) - z(1, j))^2))/(r(1, i) + r(1, j));
                if overlap(1, j) > 0
                    q = q + 1;
                    xn = xn + x(1, j) + (x(1, i) - x(1, j))*(r(1, i) + r(1, j))/sqrt((x(1, i) - x(1, j))^2 + (y(1, i) - y(1, j))^2 + (z(1, i) - z(1, j))^2);
                    yn = yn + y(1, j) + (y(1, i) - y(1, j))*(r(1, i) + r(1, j))/sqrt((x(1, i) - x(1, j))^2 + (y(1, i) - y(1, j))^2 + (z(1, i) - z(1, j))^2);
                    zn = zn + z(1, j) + (z(1, i) - z(1, j))*(r(1, i) + r(1, j))/sqrt((x(1, i) - x(1, j))^2 + (y(1, i) - y(1, j))^2 + (z(1, i) - z(1, j))^2);
                end
            end
        end
        if q >= 1
            x(1, i) = xn/q;
            y(1, i) = yn/q;
            z(1, i) = zn/q;
            
            %Dont accept if lies outside the boundary
            if (x(1, i) >= (x_container - r(1, i)))
                x(1, i) = (x_container - r(1, i));
            end
            if (y(1, i) >= (y_container - r(1, i)))
                y(1, i) = (y_container - r(1, i));
            end
            if (z(1, i) >= (z_container - r(1, i)))
                z(1, i) = (z_container - r(1, i));
            end
            if (x(1, i) <= r(1, i))
                x(1, i) = r(1, i);
            end
            if (y(1, i) <= r(1, i))
                y(1, i) = r(1, i);
            end
            if (z(1, i) <= r(1, i))
                z(1, i) = r(1, i);
            end
        end
    end
end
    
% Plotting all the particles
[X,Y,Z] = sphere(15);
for j=1:num_particles
    surf(X*r(1, j) + x(1, j), Y*r(1, j) + y(1, j), Z*r(1, j) + z(1, j));
    hold on
end

%% Finding the maximum overlap of all particles (just to check)
% First create the array of all particles
% q = 1;
% for cell = 1:total_cell_count
%     for par = 1:num_particles
%         x_total(q, 1) = x(cell, par);
%         y_total(q, 1) = y(cell, par);
%         z_total(q, 1) = z(cell, par);
%         r_total(q, 1) = r(1, par);
%         q = q + 1;
%     end
% end
% Find the maximum overlap of all particles
% maximum_overlap = 0;
% for i = 2:totall_particle_count
%     for j = 1:(i-1)
%         if (sqrt((x_total(i, 1) - x_total(j, 1))^2 + (y_total(i, 1) - y_total(j, 1))^2 + (z_total(i, 1) - z_total(j, 1))^2) < (r_total(i, 1) + r_total(j, 1)))
%             overlap(i, j) = ((r_total(i, 1) + r_total(j, 1)) - sqrt((x_total(i, 1) - x_total(j, 1))^2 + (y_total(i, 1) - y_total(j, 1))^2 + (z_total(i, 1) - z_total(j, 1))^2))/(r_total(i, 1) + r_total(j, 1));
%         end
%     end
% end
% maximum_overlap = max(overlap);

%% Finding the neighbors of each particle in the bed
% for i = 1:totall_particle_count
%     q = 1;
%     for j = 1:totall_particle_count
%         if (i ~= j) && ((sqrt((x_total(i, 1) - x_total(j, 1))^2 + (y_total(i, 1) - y_total(j, 1))^2 + (z_total(i, 1) - z_total(j, 1))^2) - (r_total(i, 1) + r_total(j, 1))) <= 0.000001)
%             neighbors(i, q) = j;
%             q = q + 1;
%         end
%     end
% end

%% Finding the number of neighbors each particle has
% sizen = size(neighbors);
% number_neighbor = zeros(sizen(2), 1);
% for i = 1:totall_particle_count
%     for j = 1:sizen(2)
%         if neighbors(i, j) == 0
%             number_neighbor(j, 1) = number_neighbor(j, 1) + 1;
%             break;
%         end
%     end
% end

%% Writing particle locations in txt files
% dlmwrite('particle_x.txt', x);
% dlmwrite('particle_y.txt', y);
% dlmwrite('particle_z.txt', z);
% dlmwrite('particle_r.txt', r);
% dlmwrite('Neighbors.txt', neighbors);

toc