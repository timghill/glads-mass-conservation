function make_trimesh()
% Script to make triangular mesh
% make the box and write it to a file which the meshing library (Triangle) knows
boundary_xy = [0, 0;
               100e3, 0;
               100e3, 25e3;
               0, 25e3];

% boundary marks>0 on edge:
bmark = [1;2;2;1];          % just a mark which is given to the nodes on the boundary
bmark_edge = [2;2;2;1];     % just a mark which is given to the edges on the boundary

% maxarea = 0.5e6; % Makes mesh with 4156 nodes
maxarea = 1e6;

% cell array holding all the meshes
mesh = make_mesh(boundary_xy, bmark, bmark_edge, maxarea);

save('mesh.mat', '-struct', 'mesh');