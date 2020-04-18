% 图像放大函数 myimresize
% 输入参数为输入矩阵，插值方法以及输出图像的大小(x为图像的width，y为图像的length)
% 输出图像大小可以是倍数也可以是放大后的尺寸
% 返回值是输出矩阵，uint8
% method 有三个可选项：'nearest'、'bilinear'、'bicubic'
function R=myimresize(I,method,x,y)
% 输入图像的列、行
[j,k,~]=size(I);
% 参数输入个数判断
if nargin==3
    x_new=floor(j*x);
    y_new=floor(k*x);
else 
    x_new=x;
    y_new=y;
end
I_dim=ndims(I);
% 最近邻算法
if method == "nearest"
    if I_dim==2
    R=nn(I,x_new,y_new);
    else 
       R1=nn(I(:,:,1),x_new,y_new);
       R2=nn(I(:,:,2),x_new,y_new);
       R3=nn(I(:,:,3),x_new,y_new);
       R= cat(3, R1, R2, R3);
    end
end
% 双线性
if method == "bilinear"
    if I_dim==2
    R=bili2(I,x_new,y_new);
    else 
       R1=bili2(I(:,:,1),x_new,y_new);
       R2=bili2(I(:,:,2),x_new,y_new);
       R3=bili2(I(:,:,3),x_new,y_new);
       R= cat(3, R1, R2, R3);
    end
end  
% 双立方
if method=="bicubic"
    if I_dim==2
    R=Bicubic(I,x_new,y_new);
    else 
       R1=Bicubic(I(:,:,1),x_new,y_new);
       R2=Bicubic(I(:,:,2),x_new,y_new);
       R3=Bicubic(I(:,:,3),x_new,y_new);
       R= cat(3, R1, R2, R3);
    end
end   
R=uint8(R);
end

% 最近邻函数
function  M=nn(I,x_new,y_new)
% 输入图像的行列
[j,k]=size(I);
% 缩放因子
x_scale = x_new./(j);
y_scale = y_new./(k);
% 输出矩阵
M = zeros(x_new,y_new);
% 核心代码
for cou1 = 1:x_new
    for cou2 = 1:y_new
        % 取最近的点，实际上用ceil更合适，而非round,同时解决了边界范围问题
        % 映射点落在哪个像素就用哪个像素，也符合最近邻的思想
        M(cou1,cou2) = I(ceil(cou1./x_scale),ceil(cou2./y_scale));
    end
end
end

% 双线性
function M=bili2(I,x_new,y_new)
% 输入图像的行列
[j,k]=size(I);
% 缩放因子
x_scale = x_new./(j-1);
y_scale = y_new./(k-1);
% 输出矩阵
M = zeros(x_new,y_new);
% 核心代码
for cou1 = 0:x_new-1
    for cou2 = 0:y_new-1
        % 输出图像映射到输入图像的亚像素位置
        o1=cou1./x_scale;
        o2=cou2./y_scale;
        % 计算因子
        W = 1-o1+floor(o1);
        H = 1-o2+floor(o2);
        % 周围4个点
        I11 = I(1+floor(o1),1+floor(o2));
        I12 = I(1+ceil(o1),1+floor(o2));
        I21 = I(1+floor(o1),1+ceil(o2));
        I22 = I(1+ceil(o1),1+ceil(o2));
        % 计算公式
        M(cou1+1,cou2+1) = (1-W).*(1-H).*I22 + (W).*(1-H).*I21 + (1-W).*(H).*I12 + (W).*(H).*I11;
    end
end
end

% 双立方插值
function M=Bicubic(I,x_new,y_new)
% 输入图像的行列
[j,k]=size(I);
% 缩放因子
s1=x_new/j;
s2=y_new/k;
% 输入图像的padding矩阵：j+3行，k+3列
I_pad=[I(1,1),I(1,:),repmat(I(1,end),1,2);...
       I(:,1),I,repmat(I(:,end),1,2);...
       repmat(I(end,1),2,1),repmat(I(end,:),2,1),repmat(I(end,end),2,2)];
% uint8转换为double
I_pad=double(I_pad);
% 多项式p计算函数,I是列向量
N=[0 2 0  0;...
    -1 0 1 0;...
    2 -5 4 -1;...
    -1 3 -3 1 ];
p=@(d,I) 0.5*[1 d d^2 d^3]*N*I;
% 输入图像的亚像素位置，映射到输出图像
rf=1:x_new;
cf=1:y_new;
rf_to_I=rf./s1;
cf_to_I=cf./s2;
% r0,c0：像素的行列索引
r0=floor(rf_to_I);
c0=floor(cf_to_I);
% 边界范围限制
r0=(r0<1).*1+(r0>=1).*r0;
c0=(c0<1).*1+(c0>=1).*c0;
% 多项式p中的定位元素
dr=rf_to_I-r0;
dc=cf_to_I-c0;
% K为列向量
K=zeros(4,1);
% 输出图像矩阵
M=zeros(x_new,y_new);
for row=1:x_new
    for column=1:y_new
        % 先沿列计算
        for m=0:3
            K(m+1)=p(dr(row),I_pad((r0(row):r0(row)+3),c0(column)+m)); 
        end
        % 再按照行计算输出像素值
        M(row,column)=p(dc(column),K);
    end
end
end
