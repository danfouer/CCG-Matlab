% 这是一个用于将网格数据写入OFF格式文件的MATLAB函数。OFF是一种常见的三维网格文件格式，包含了网格的拓扑结构
% 和几何信息。该函数可以将面片连接信息和顶点位置信息写入OFF文件中，并可选地写入顶点或面片的颜色信息。
% 
% 函数的输入参数为OFF文件名、面片连接信息、顶点位置信息和颜色信息（可选），输出参数为空。函数首先打开文件，
% 写入文件头信息，然后将顶点位置信息和面片连接信息写入文件中。在写入顶点位置信息时，函数使用了dlmwrite函数，
% 将每个顶点的坐标写入文件中。在写入面片连接信息时，函数同样使用了dlmwrite函数，将每个面片的连接信息写入文件中。
% 如果颜色信息不为空，则将颜色信息与顶点位置信息或面片连接信息合并后一起写入文件中。
% 
% 该函数的作者是 Meng Bin，版权归香港中文大学数学系计算几何组所有。
%% write off 
% Write mesh data to OFF format mesh file
%  
%% Syntax
%   write_off(filename,face,vertex,color)
%   write_off(filename,face,vertex)
%
%% Description
%  filename: string, file to read.
%  face    : double array, nf x 3 array specifying the connectivity of the mesh.
%  vertex  : double array, nv x 3 array specifying the position of the vertices.
%  color   : double array, nv x 3 or nf x 3 array specifying the color of the vertices or faces.
%
%%  Example
%   write_off('temp.off',face,vertex);
%   write_off('temp.off',face,vertex,clor);
%
%% Contribution
%  Author : Meng Bin
%  Created: 2014/03/05
%  Revised: 2014/03/07 by Meng Bin, Block write to enhance writing speed.
%  Revised: 2014/03/17 by Meng Bin, modify doc format
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui

function write_off(filename,face,vertex,color)

if nargin < 4
    color = [];
end

fid = fopen(filename,'wt');
if( fid==-1 )
    error('Can''t open the file.');
end

nvert = size(vertex, 1);
nface = size(face, 1);
nvert_face = size(face, 2);

ncolor =0;
if ~isempty(color)
    ncolor = size(color, 1);
end

fprintf (fid, 'OFF\n');
fprintf (fid, '%d %d %d\n',nvert, nface, 0);

if nvert == ncolor 
    vertex = [vertex';color']';
end
if nface == ncolor && nvert ~= ncolor
    face =[zeros(1,nface)+nvert_face; face'-1;color']';
else
    face =[zeros(1,nface)+nvert_face;face'-1]';
end
dlmwrite(filename,vertex,'-append',...
         'delimiter',' ',...
         'precision', 6,...
         'newline','pc');

dlmwrite(filename,face,'-append',...
         'delimiter',' ',...
         'newline','pc');

fclose(fid);
