% 这段 MATLAB 代码实现了计算三角网格的边权重。具体实现过程如下：
% 
% 首先，从输入的 mesh 结构中获取 face、edge 和 vert 矩阵，分别表示三角网格的面、边和顶点。
% 同时，获取 eif 矩阵，它是一个大小为 ne x 2 的矩阵，表示每条边所在的两个面的编号。如果一条边只属于一个面，
% 则另一个面的编号为 0。
% 
% 然后，根据 eif 矩阵将边分为两类：属于两个面的内部边和只属于一个面的边。对于内部边，
% 计算它所在的两个面的第三个顶点 ev1 和 ev2，以及它所对应的两个三角形的内角余切值 ct1 和 ct2。
% 将 ct1 和 ct2 分别加到对应边的权重 ew 中。对于只属于一个面的边，不需要计算权重。
% 
% 接着，将属于两个面的内部边的权重除以 2，因为每条边会被计算两次。
% 
% 最后，将计算得到的边权重 ew 存储在 mesh 结构中，并将其作为函数的输出。此外，函数还定义了一个内部函数 cot2，
% 用于计算三角形的内角余切值。
function ew = edge_weight(mesh)
% compute edge weight
face = mesh.face;
edge = mesh.edge;
vert = mesh.vert;
% eif indicate which face this edge belongs to, one row for one edge, -1
% means boundary edge
eif = mesh.eif;

ne = size(edge,1);
ew = zeros(ne,1);

% ev1 is the third vert in a triangle: v3 = (v1+v2+v3) - (v1+v2)
ind1 = eif(:,1)>0;
ev1 = sum(face(eif(ind1,1),:),2) - sum(edge(ind1,:),2);
ct1 = cot2(vert(ev1,:),vert(edge(ind1,1),:),vert(edge(ind1,2),:));
ew(ind1) = ew(ind1) + ct1;
% ev2 is similar to ev1
ind2 = eif(:,2)>0;
ev2 = sum(face(eif(ind2,2),:),2) - sum(edge(ind2,:),2);
ct2 = cot2(vert(ev2,:),vert(edge(ind2,1),:),vert(edge(ind2,2),:));
ew(ind2) = ew(ind2) + ct2;
ew(ind1&ind2) = ew(ind1&ind2)/2;
% store ew for later use
mesh.ew = ew;

function ct = cot2(pi,pj,pk)
a = sqrt(dot(pj-pk,pj-pk,2));
b = sqrt(dot(pk-pi,pk-pi,2));
c = sqrt(dot(pi-pj,pi-pj,2));
cs = (b.*b+c.*c-a.*a)./(2.*b.*c);
ss2 = 1-cs.*cs;
ss2(ss2<0) = 0;
ss2(ss2>1) = 1;
ss = sqrt(ss2);
ct = cs./ss;