function [load_x,load_y] = ComputeLoadGraphics(y,l_rope)

x_com   = y(1);
y_com   = y(3);
x_load  = y(15);
y_load  = y(17);

l_load = 1.0;
h_load = 1.0;

l_text = 0.08;
% Define vertices of box 
load_x = [x_load-0.5*l_load,    x_load+0.5*l_load... 
         x_load+0.5*l_load,    x_load-0.5*l_load];
load_y = [y_load-0.5*h_load,    y_load-0.5*h_load...
        y_load+0.5*h_load,     y_load+0.5*h_load];

d_quad_load  = sqrt( (y(1)-y(15))^2 + (y(3)-y(17))^2);

end

