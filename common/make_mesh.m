% 这是一个用于构造网格的MATLAB函数。该函数的输入参数为一个三角面片的连接关系face和顶点坐标vert，
% 输出参数为一个包含网格信息的结构体mesh。
% 
% 函数的具体实现如下：
% 
% 1. 根据输入的面片连接关系face，使用edges函数计算出网格的边的连接关系edge和每条边所属的两个面片的编号eif。
% 
% 2. 根据输入的面片连接关系face，使用halfedges函数计算出网格的半边连接关系he和每条半边所属的面片的编号heif。
% 
% 3. 计算网格的面片数nf、顶点数nv、边数ne和半边数nh，并将这些信息存储在mesh结构体中。
% 
% 4. 将计算得到的面片连接关系face、顶点坐标vert、边的连接关系edge、每条边所属的两个面片的编号eif、
% 半边连接关系he和每条半边所属的面片的编号heif存储在mesh结构体中。
% 
% 5. 计算网格的边界点，并将其存储在mesh结构体的bd字段中。其中，bde表示所有在边界上的边，
% unique(bde(:))表示所有在边界上的点。
function mesh = make_mesh(face, vert)

[edge,eif] = edges(face);
[he,heif] = halfedges(face);

nf = size(face,1);
nv = size(vert,1);
ne = size(edge,1);
nh = size(he,1);

mesh.nf = nf;
mesh.nv = nv;
mesh.ne = ne;
mesh.nh = nh;

mesh.face = face;
mesh.vert = vert;
mesh.edge = edge;
mesh.eif = eif;
mesh.halfedge = he;
mesh.heif = heif;

% bd is unordered, mainly used to test if vert is on boundary
ind = eif(:,1)>0 & eif(:,2)>0;
bde = edge(~ind,:);
mesh.bd = unique(bde(:));
