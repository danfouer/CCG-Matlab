% 这是一个计算网格上最短路径的MATLAB函数。该函数的输入参数为一个包含网格信息的结构体mesh
% 和一个包含路径上顶点编号的向量cc，输出参数为一个包含路径上顶点编号的向量path。
% 
% 函数的具体实现如下：
% 
% 1. 从mesh结构体中获取顶点坐标vert和边的连接关系edge，并计算出每条边的长度el。
% 
% 2. 使用graph函数构造一个无向图G，其中顶点为所有的顶点，边为所有的边，边的权重为对应的边的长度el。
% 
% 3. 对于路径上相邻的两个顶点cc(i-1)和cc(i)，使用shortestpath函数计算它们之间的最短路径pathi，
% 并将其存储在pathi中。
% 
% 4. 将所有的pathi拼接起来，得到最终的路径path。
% 
% 5. 将路径path转化为列向量，并将其作为函数的输出。
function path = shortest_path(mesh,cc)
vert = mesh.vert;
edge = mesh.edge;
de = vert(edge(:,1),:)-vert(edge(:,2),:);
el = sqrt(dot(de,de,2));
G = graph([edge(:,1);edge(:,2)],[edge(:,2);edge(:,1)],[el;el]);
path = [];
for i = 2:length(cc)
    pathi = shortestpath(G,cc(i-1),cc(i));
    path = [path,pathi(1:end-1)];
end
path = [path,cc(end)];
path = path(:);
