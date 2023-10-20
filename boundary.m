function [u_new, v_new] = boundary(x, y, W, H, u, v, V) 
%Function will:
%find distance between given point (individual) and the boundaries
%If individual less than 2m from a boundary (wall), a new direction (theta) is randomly generated
%Direction will be generated in such a way, that will return individual back to domain
%Finally, velocity components will be updated and outputed 

%Find distance between given point and the boundaries
b_l = x; %left boundary
b_r = W - x; %right boundary
b_t = H - y; %top boundary
b_b = y; %bottom boundary

%Set output velocities and position equal to input, to make sure that
%function outputs variables if no modification of variables is required
u_new = u;
v_new = v;

%If individual less than 2m from a wall, generate new direction and update velocity components 
if b_l < 2
    theta = pi*rand-pi/2; %New random direction between -90° and 90°
    u_new = V*cos(theta);
    v_new = V*sin(theta);
end

if b_r < 2
    theta = pi*rand+pi/2; %New random direction between 90° and 270°
    u_new = V*cos(theta);
    v_new = V*sin(theta);
end

if b_t < 2
    theta = pi*rand+pi; %New random direction between 180° and 360°
    v_new = V*sin(theta);
    u_new = V*cos(theta);
end

if b_b < 2
    theta = pi*rand; %New random direction between 0° and 180°
    v_new = V*sin(theta);
    u_new = V*cos(theta);
end
end