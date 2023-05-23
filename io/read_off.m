% 这是一个用于读取OFF格式文件的MATLAB函数。该函数支持读取网格的连接关系和顶点位置，其他数据将被丢弃。
% 
% 函数的语法为：
% 
% [face,vertex] = read_off(filename)
% [face,vertex,extra] = read_off(filename)
% 
% 其中，`filename` 是要读取的文件名，`face` 是一个大小为 `nf x 3` 的双精度数组，表示网格的连接关系，`vertex` 是一个大小为 `nv x 3` 的双精度数组，表示顶点的位置，`extra` 是一个结构体，包含除了 `face` 和 `vertex` 以外的所有数据。
% 
% 函数的实现过程如下：
% 
% 1. 打开文件，判断是否为 OFF 文件。
% 2. 读取顶点和面的数量。
% 3. 读取顶点数据。
% 4. 读取面数据。
% 5. 如果顶点数据包含颜色信息，则将颜色信息存入 `extra.Vertex_color` 中。
% 6. 将面数据的索引加 1，因为 OFF 文件中的索引从 0 开始。
% 7. 如果面数据包含颜色信息，则将颜色信息存入 `extra.Face_color` 中。
% 8. 返回 `face`、`vertex` 和 `extra`。
% 
% 需要注意的是，该函数假设 OFF 文件中的面都是三角形，因此只返回三角形网格的连接关系。
%% read_off
% Read mesh data from OFF file.
% 
%% Syntax
%   [face,vertex] = read_off(filename)
%   [face,vertex,extra] = read_off(filename)
%
%% Description
%  filename: string, file to read.
%
%  face  : double array, nf x 3, connectivity of the mesh.
%  vertex: double array, nv x 3, position of the vertices.
%  extra : struct, anything other than face and vertex are included.
%
%% Example
%   [face,vertex] = read_off('torus.off');
%   [face,vertex,extra] = read_off('torus.off');
%
%% Contribution
%  Author : Wen Cheng Feng
%  Created: 2014/03/27
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui

function [face,vertex,extra] = read_off(filename)
fid = fopen(filename,'r');
if( fid == -1 )
    error('Can''t open the file.');
end

% determine if this is an OFF file
line = get_next_line(fid);
type = sscanf(line,'%s');
if ~strcmp(type,'OFF') && ~strcmp(type,'COFF') && ~strcmp(type,'CNOFF')
    fclose(fid);
    error('Not a valid OFF file.');    
end
% read vertex and face number
line = get_next_line(fid);
nvf = sscanf(line,'%d');
nv = nvf(1);
nf = nvf(2);

% read vertex
k = 0;
vertex = [];
line = get_next_line(fid);
[format,sz] = get_format(line);
fseek(fid,-length(line),'cof');
while k < nv
    A = fscanf(fid,format,(nv-k)*sz);
    vertex = [vertex;reshape(A,sz,length(A)/sz)'];
    k = k + length(A)/sz;
    line = get_next_line(fid);
    if line == -1 % end of file
        if k ~= nv
            error('vertex data is not correct.');
        end
    end
    fseek(fid,-length(line),'cof');
end

% read face
k = 0;
face = [];
line = get_next_line(fid);
[format,sz] = get_format(line);
fseek(fid,-length(line),'cof');
while k < nf
    A = fscanf(fid,format,(nf-k)*sz);
    face = [face;reshape(A,sz,length(A)/sz)'];
    k = k + length(A)/sz;
    line = get_next_line(fid);
    if line == -1 % end of file
        if k ~= nf
            error('face data is not correct.');
        end
    end
    fseek(fid,-length(line),'cof');
end

extra = [];
if size(vertex,2) > 3
    extra.Vertex_color = vertex(:,4:end);
    vertex = vertex(:,1:3);
end
face(:,2:4) = face(:,2:4) + 1;
if size(face,2) > 4
    extra.Face_color = face(:,5:end);    
end
face = face(:,2:4);

fclose(fid);

function line = get_next_line(fid)
% read next line, skip comment, blank line and eof
line = fgets(fid);
if line == -1    
    return
end
while isempty(strtrim(line)) || line(1) == '#'
    line = fgets(fid);
    if line == -1
        return;
    end
end

function [format,sz] = get_format(str)
% determine the format of input str
sn = '[\-+]?(?:\d*\.|)\d+\.?(?:[eE][\-+]?\d+|)'; % match number
format = regexprep(str,sn,'%f');
[~,splitstr] = regexp(str,sn,'match');
sz = length(splitstr);