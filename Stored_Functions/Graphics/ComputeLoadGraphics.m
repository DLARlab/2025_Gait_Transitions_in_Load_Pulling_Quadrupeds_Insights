function [load_x,load_y] = ComputeLoadGraphics(y)


x_load  = y(15);
y_load  = y(17);

l_load = y_load;
h_load = y_load;


% Define vertices of box 
load_x = [x_load-l_load,    x_load+l_load... 
         x_load+l_load,    x_load-l_load];
load_y = [y_load-h_load,    y_load-h_load...
        y_load+h_load,     y_load+h_load];

end

