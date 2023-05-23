% 该函数用于计算调和一形式。
% 
% 输入参数：
% - mesh: 包含网格信息的结构体。
% - w: double类型，n x 3，表示一个闭合的一形式。
% 
% 输出参数：
% - dh: double类型，n x 3，表示在w的同调类中的调和一形式。
% 
% 函数的实现过程：
% 1. 使用exterior_co_derivative函数计算w的外共轭。
% 2. 使用laplace_beltrami函数计算拉普拉斯-贝尔特拉米算子L。
% 3. 将L的第一行第一列元素减1，以满足调和一形式的定义。
% 4. 使用反斜杠运算符求解方程Lh = -delta_w，得到调和函数h。
% 5. 使用exterior_derivative函数计算dh，它是一个闭合的调和一形式。
% 6. 将w和dh相加，得到在w的同调类中的调和一形式。
function dh = harmonic_form(mesh, w)
% given a closed one form w, compute the harmonic one form in the
% homological class of w
% find function h such that delta(w + dh) = 0, delta(dh) = -delta_w
% since d(w + dh) = dw + d(dh) = 0, so w+dh is harmonic form
delta_w = exterior_co_derivative(mesh,w,1);
L = laplace_beltrami(mesh);
L(1,1) = L(1,1)-1;
h = -L\delta_w;
dh = exterior_derivative(mesh,h,0);
% dh is closed harmonic one form
dh = w + dh;
