function output = C_LGCV_to_ECEF(GS_ECEF)
%test
% clear
% clc
% lgcv = [-379.21;-348.80;-865.13];
% LLH = [-33, 151, 52];
%result = [-4678800; 2593900; -3474600];
LLH = ECEF_to_LLH(GS_ECEF(1),GS_ECEF(2),GS_ECEF(3));
%latitude
lambda_i = LLH(1) * pi/180;
%longitude
phi = LLH(2) * pi/180;
% %height with earth radius
% R = LLH(3);
Clgcv2ecef = [-sin(lambda_i)*cos(phi), -sin(phi), -cos(lambda_i)*cos(phi);...
               -sin(lambda_i)*sin(phi),  cos(phi), -cos(lambda_i)*sin(phi);...
                cos(lambda_i)         ,  0       , -sin(lambda_i)];
output = Clgcv2ecef;
end