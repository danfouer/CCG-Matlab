% 这是一个用于将网格数据写入OBJ格式文件的MATLAB函数。OBJ是一种常见的三维网格文件格式，
% 包含了网格的拓扑结构和几何信息。该函数可以将面片连接信息和顶点位置信息写入OBJ文件中。
% 
% 函数的输入参数为OBJ文件名、面片连接信息和顶点位置信息，输出参数为空。函数首先打开文件，
% 写入文件头信息，然后将顶点位置信息和面片连接信息写入文件中。在写入顶点位置信息时，
% 函数使用了fprintf函数，将每个顶点的坐标写入文件中。在写入面片连接信息时，函数同样使用了fprintf函数，
% 将每个面片的连接信息写入文件中。
% 
% 该函数的作者是 Meng Bin，版权归香港中文大学数学系计算几何组所有。
%% write obj 
% Write mesh data to OBJ format mesh file
%
%% Syntax
%   write_obj(filename,face,vertex)
%
%% Description
%  filename: string, file to read.
%  face    : double array, nf x 3 array specifying the connectivity of the mesh.
%  vertex  : double array, nv x 3 array specifying the position of the vertices.
%  color   : double array, nv x 3 or nf x 3 array specifying the color of the vertices or faces.
%
%% Example
%   write_obj('cube.obj',face,vertex);
%
%% Contribution
%  Author : Meng Bin
%  History: 2014/03/05 file created
%  Revised: 2014/03/07 by Meng Bin, Block write to enhance writing speed
%  Revised: 2014/03/17 by Meng Bin, modify doc format
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui

function write_obj(filename,face,vertex)

fid = fopen(filename,'w');
if( fid==-1 )
    error('Can''t open the file.');
end

% write logo
fprintf (fid, '#Generated by geometric processing package.\n');

% write vertex
fprintf (fid, 'v %.6f %.6f %.6f\n',vertex');

% write face
fprintf (fid, 'f %d %d %d\n',face');

fclose(fid);
