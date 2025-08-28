% The code can be used to simulate the velocity and strain rate values by using only the geometric dimensions and the coefficient of friction.
% Operates on MATLAB R2023a and GNU OCTAVE 9.2
% The model is published in the following article:
% Bátorfi, J. G., & Sidor, J. J. (2025). Computationally Effective Modeling of Cold Rolling: Application to Al Alloys. Metals, 15(1), 11.
% https://doi.org/10.3390/met15010011


% Clear workspace
clear;
clc;
% INPUTS

% Geometrical parameters
h_1=2.26;     % half-thickness after rolling (mm)
h_0=2.51;     % half-thickness before rolling (mm)
R=75;         % radius of the roll (mm)
omega=1.101;  % angular velocoty of rolls (rad/s)
mu=0.15;      % friction coefficient
n=100;        % flowline exponent

% Numerical parameters for plot and calculation
N=200;        % number of points along to the flowline
M=10;         % number of paralel flowlines in half of the thickness

% CALCULATIONS, DO NOT MODIFY THE CALCULATION BELOW THIS LINE

% Initial calculations
d=R*sin(acos((R+h_1-h_0)/R));
L_d=d;
h_i=2*h_0;
h_f=2*h_1;
mu_min=0.5*sqrt(h_f/R)*(log(h_i/h_f)+0.25*sqrt((h_i-h_f)/R))/(atan(h_i/h_f-1));
mu_fact=mu/mu_min;
a_21_pre=0.1626213*(h_0/R)^1.688062*(h_1/R)^(-0.131232)*((h_0-h_1)/R)^(-0.97944);
a_21_exp=43.56109*(h_0/R)^(-3.31412)*(h_1/R)^2.509572*((h_0-h_1)/R)^0.986799;
a_21=a_21_pre*(mu/mu_min)^a_21_exp;
a_22=0.986307*(h_0/R)^(-0.1132);
a_23_pre=0.308969*(h_0/R)^(-4.21639)*(h_1/R)^4.132604*((h_0-h_1)/R)^0.323943;
a_23_exp=0.0373779*(h_0/R)^6.493182*(h_1/R)^(-4.173077)*((h_0-h_1)/R)^(-2.16706);
a_23=a_23_pre*(mu/mu_min )^a_23_exp;
a_24=3.025;
a_20=-(a_21-a_23);
H_0=2*sqrt(R/h_f)*atan(sqrt((h_i-h_f)/h_f));
H_n=0.5*(H_0-1/mu*log(h_i/h_f));
fi_n=sqrt(h_f/R)*tan(sqrt(h_f/R)*H_n/2);
L_dN=L_d-R*sin(fi_n);
q_1=1.042896+0.308671*(d/R)^0.01738*eps^(-0.00419);
q_2=12.3624-0.02363*(d/R)^(-0.18077)*eps^0.170827;
q_3=0.148264+0.000927*(d/R)^0.758336*eps^(-0.19655);
q_4=0.322031+0.086773*(d/R)^2.109206*eps^0.51112;
l_d=L_d;
h_n=(h_0-h_1)/(q_1*exp(q_2*(L_dN/d-q_3))+1)^q_4+h_1;
h_dot_n=-((h_0-h_1)*q_1*q_2*q_4*exp(q_2*(L_dN/d-q_3))*(q_1*exp(q_2*(L_dN/d-q_3))+1)^((-q_4)-1))/d;
a_12=-1/h_dot_n;
f_0=1/h_n*h_1;
v_0=R*omega;
lam=v_0*f_0*h_0;
v_out=lam/h_1;

% Initializing vectors and matrices
z_n_ser=0.001:0.999/M:1;
x_ser=-0.5*d:1.6*d/N:1.5*d;
z=zeros(1.25*N,M+1);
x=zeros(1.25*N,M+1);
z_s=zeros(1.25*N,M+1);
z_n=zeros(1.25*N,M+1);
v_x=zeros(1.25*N,M+1);
v_z=zeros(1.25*N,M+1);
v_tot=zeros(1.25*N,M+1);
L_xx=zeros(1.25*N,M+1);
L_xz=zeros(1.25*N,M+1);
L_zx=zeros(1.25*N,M+1);
L_zz=zeros(1.25*N,M+1);
D_xx=zeros(1.25*N,M+1);
D_xz=zeros(1.25*N,M+1);
D_zx=zeros(1.25*N,M+1);
D_zz=zeros(1.25*N,M+1);
D_xx_real=zeros(1.25*N,M+1);
D_xz_real=zeros(1.25*N,M+1);
D_zx_real=zeros(1.25*N,M+1);
D_zz_real=zeros(1.25*N,M+1);
t=zeros(1.25*N,M+1);
dx=zeros(1.25*N,M+1);
h=zeros(1.25*N,M+1);
h_dot_x=zeros(1.25*N,M+1);
h_dot_z=zeros(1.25*N,M+1);
h_ddot_xx=zeros(1.25*N,M+1);
h_ddot_xz=zeros(1.25*N,M+1);
h_ddot_zz=zeros(1.25*N,M+1);
h_dddot_xxx=zeros(1.25*N,M+1);
h_dddot_xxz=zeros(1.25*N,M+1);
h_dddot_xzz=zeros(1.25*N,M+1);
h_dddot_zzz=zeros(1.25*N,M+1);
phi=zeros(1.25*N,M+1);
eps_dot_11=zeros(1.25*N,M+1);
eps_dot_13=zeros(1.25*N,M+1);

% Calculation of coordinates and properties for the points of flowlines
for i=1:(1.25*N)
  for k=1:(M+1)

% Calculation of coordinates for further calculation
    x(i,k)=x_ser(i);
    z_n(i,k)=z_n_ser(k);
    h(i,k)=(h_0-h_1)/(q_1*exp(q_2*(x(i,k)/d-q_3))+1)^q_4+h_1;
    z(i,k)=h(i,k)*z_n(i,k);

% Calculation of parameters
    h_dot_x(i,k)=-((h_0-h_1)*q_1*q_2*q_4*exp(q_2*(x(i,k)/d-q_3))*(q_1*exp(q_2*(x(i,k)/d-q_3))+1)^((-q_4)-1))/d;
    h_dot_z(i,k)=0;
    h_ddot_xx(i,k)=(-((h_0-h_1)*q_1^2*q_2^2*((-q_4)-1)*q_4*(q_1*exp(q_2*(x(i,k)/d-q_3))+1)^((-q_4)-2)*exp(2*q_2*(x(i,k)/d-q_3)))/d^2)-((h_0-h_1)*q_1*q_2^2*q_4*exp(q_2*(x(i,k)/d-q_3))*(q_1*exp(q_2*(x(i,k)/d-q_3))+1)^((-q_4)-1))/d^2;;
    h_ddot_zz(i,k)=0;
    h_ddot_xz(i,k)=0;
    h_dddot_xxx(i,k)=-((h_0-h_1)*q_1^3*q_2^3*(-q_4-2)*(-q_4-1)*q_4*(q_1*exp(q_2*(x(i,k)/d-q_3))+1)^(-q_4-3)*exp(3*q_2*(x(i,k)/d-q_3)))/d^3-(3*(h_0-h_1)*q_1^2*q_2^3*(-q_4-1)*q_4*(q_1*exp(q_2*(x(i,k)/d-q_3))+1)^(-q_4-2)*exp(2*q_2*(x(i,k)/d-q_3)))/d^3-((h_0-h_1)*q_1*q_2^3*q_4*exp(q_2*(x(i,k)/d-q_3))*(q_1*exp(q_2*(x(i,k)/d-q_3))+1)^(-q_4-1))/d^3;
    h_dddot_xxz(i,k)=0;
    h_dddot_xzz(i,k)=0;
    h_dddot_zzz(i,k)=0;
    z_s(i,k)=(z(i,k)*((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+1))/h(i,k);
    phi(i,k)=z_s(i,k);

% Calculation of velocities
    v_x(i,k)=lam*((z(i,k)*((2*a_12*h_dot_x(i,k)*h_ddot_xz(i,k)+h_ddot_xz(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+(a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((a_21*a_22*h(i,k)*(z(i,k)/h(i,k))^a_22*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k)-(a_23*a_24*h(i,k)*(z(i,k)/h(i,k))^a_24*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k))))/h(i,k)-(h_dot_z(i,k)*z(i,k)*((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+1))/h(i,k)^2+((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+1)/h(i,k));
    v_z(i,k)=-lam*((z(i,k)*((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((a_23*a_24*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_24)/h(i,k)-(a_21*a_22*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_22)/h(i,k))+(2*a_12*h_dot_x(i,k)*h_ddot_xx(i,k)+h_ddot_xx(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)))/h(i,k)-(h_dot_x(i,k)*z(i,k)*((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+1))/h(i,k)^2);
    v_tot(i,k)=sqrt(v_x(i,k)^2+  v_z(i,k)^2);

   % Calculation of strain velocities
        L_xx(i,k)=lam*((z(i,k)*((2*a_12*h_dot_x(i,k)*h_ddot_xz(i,k)+h_ddot_xz(i,k))*((a_23*a_24*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_24)/h(i,k)-(a_21*a_22*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_22)/h(i,k))+(2*a_12*h_dot_x(i,k)*h_dddot_xxz(i,k)+h_dddot_xxz(i,k)+2*a_12*h_ddot_xz(i,k)*h_ddot_xx(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+(a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((-(a_23*a_24*h(i,k)*(z(i,k)/h(i,k))^a_24*((2*h_dot_x(i,k)*h_dot_z(i,k)*z(i,k))/h(i,k)^3-(h_ddot_xz(i,k)*z(i,k))/h(i,k)^2-h_dot_x(i,k)/h(i,k)^2))/z(i,k))+(a_21*a_22*h(i,k)*(z(i,k)/h(i,k))^a_22*((2*h_dot_x(i,k)*h_dot_z(i,k)*z(i,k))/h(i,k)^3-(h_ddot_xz(i,k)*z(i,k))/h(i,k)^2-h_dot_x(i,k)/h(i,k)^2))/z(i,k)+(a_23*a_24^2*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_24*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k)-(a_23*a_24*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_24*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k)-(a_21*a_22^2*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_22*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k)+(a_21*a_22*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_22*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k))+(2*a_12*h_dot_x(i,k)*h_ddot_xx(i,k)+h_ddot_xx(i,k))*((a_21*a_22*h(i,k)*(z(i,k)/h(i,k))^a_22*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k)-(a_23*a_24*h(i,k)*(z(i,k)/h(i,k))^a_24*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k))))/h(i,k)-(h_dot_z(i,k)*z(i,k)*((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((a_23*a_24*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_24)/h(i,k)-(a_21*a_22*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_22)/h(i,k))+(2*a_12*h_dot_x(i,k)*h_ddot_xx(i,k)+h_ddot_xx(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)))/h(i,k)^2+((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((a_23*a_24*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_24)/h(i,k)-(a_21*a_22*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_22)/h(i,k))+(2*a_12*h_dot_x(i,k)*h_ddot_xx(i,k)+h_ddot_xx(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20))/h(i,k)-(h_dot_x(i,k)*z(i,k)*((2*a_12*h_dot_x(i,k)*h_ddot_xz(i,k)+h_ddot_xz(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+(a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((a_21*a_22*h(i,k)*(z(i,k)/h(i,k))^a_22*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k)-(a_23*a_24*h(i,k)*(z(i,k)/h(i,k))^a_24*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k))))/h(i,k)^2+(2*h_dot_x(i,k)*h_dot_z(i,k)*z(i,k)*((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+1))/h(i,k)^3-(h_ddot_xz(i,k)*z(i,k)*((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+1))/h(i,k)^2-(h_dot_x(i,k)*((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+1))/h(i,k)^2);
        L_xz(i,k)=lam*((z(i,k)*((2*a_12*h_dot_x(i,k)*h_dddot_xzz(i,k)+h_dddot_xzz(i,k)+2*a_12*(h_ddot_xz(i,k))^2)*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+(a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((-(a_23*a_24^2*h(i,k)^2*(z(i,k)/h(i,k))^a_24*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2)^2)/z(i,k)^2)+(a_21*a_22^2*h(i,k)^2*(z(i,k)/h(i,k))^a_22*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2)^2)/z(i,k)^2-(a_23*a_24*h(i,k)*(z(i,k)/h(i,k))^a_24*((-(h_ddot_zz(i,k)*z(i,k))/h(i,k)^2)+(2*(h_dot_z(i,k))^2*z(i,k))/h(i,k)^3-(2*h_dot_z(i,k))/h(i,k)^2))/z(i,k)+(a_21*a_22*h(i,k)*(z(i,k)/h(i,k))^a_22*((-(h_ddot_zz(i,k)*z(i,k))/h(i,k)^2)+(2*(h_dot_z(i,k))^2*z(i,k))/h(i,k)^3-(2*h_dot_z(i,k))/h(i,k)^2))/z(i,k)-(a_23*a_24*h_dot_z(i,k)*(z(i,k)/h(i,k))^a_24*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k)+(a_23*a_24*h(i,k)*(z(i,k)/h(i,k))^a_24*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k)^2+(a_21*a_22*h_dot_z(i,k)*(z(i,k)/h(i,k))^a_22*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k)-(a_21*a_22*h(i,k)*(z(i,k)/h(i,k))^a_22*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k)^2)+2*(2*a_12*h_dot_x(i,k)*h_ddot_xz(i,k)+h_ddot_xz(i,k))*((a_21*a_22*h(i,k)*(z(i,k)/h(i,k))^a_22*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k)-(a_23*a_24*h(i,k)*(z(i,k)/h(i,k))^a_24*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k))))/h(i,k)-(2*h_dot_z(i,k)*z(i,k)*((2*a_12*h_dot_x(i,k)*h_ddot_xz(i,k)+h_ddot_xz(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+(a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((a_21*a_22*h(i,k)*(z(i,k)/h(i,k))^a_22*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k)-(a_23*a_24*h(i,k)*(z(i,k)/h(i,k))^a_24*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k))))/h(i,k)^2+(2*((2*a_12*h_dot_x(i,k)*h_ddot_xz(i,k)+h_ddot_xz(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+(a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((a_21*a_22*h(i,k)*(z(i,k)/h(i,k))^a_22*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k)-(a_23*a_24*h(i,k)*(z(i,k)/h(i,k))^a_24*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k))))/h(i,k)-(h_ddot_zz(i,k)*z(i,k)*((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+1))/h(i,k)^2+(2*(h_dot_z(i,k))^2*z(i,k)*((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+1))/h(i,k)^3-(2*h_dot_z(i,k)*((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+1))/h(i,k)^2);
        L_zx(i,k)=-lam*((z(i,k)*((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((a_23*a_24*h_ddot_xx(i,k)*(z(i,k)/h(i,k))^a_24)/h(i,k)-(a_23*a_24^2*(h_dot_x(i,k))^2*(z(i,k)/h(i,k))^a_24)/h(i,k)^2-(a_23*a_24*(h_dot_x(i,k))^2*(z(i,k)/h(i,k))^a_24)/h(i,k)^2-(a_21*a_22*h_ddot_xx(i,k)*(z(i,k)/h(i,k))^a_22)/h(i,k)+(a_21*a_22^2*(h_dot_x(i,k))^2*(z(i,k)/h(i,k))^a_22)/h(i,k)^2+(a_21*a_22*(h_dot_x(i,k))^2*(z(i,k)/h(i,k))^a_22)/h(i,k)^2)+2*(2*a_12*h_dot_x(i,k)*h_ddot_xx(i,k)+h_ddot_xx(i,k))*((a_23*a_24*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_24)/h(i,k)-(a_21*a_22*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_22)/h(i,k))+(2*a_12*h_dot_x(i,k)*h_dddot_xxx(i,k)+h_dddot_xxx(i,k)+2*a_12*(h_ddot_xx(i,k))^2)*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)))/h(i,k)-(2*h_dot_x(i,k)*z(i,k)*((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((a_23*a_24*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_24)/h(i,k)-(a_21*a_22*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_22)/h(i,k))+(2*a_12*h_dot_x(i,k)*h_ddot_xx(i,k)+h_ddot_xx(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)))/h(i,k)^2-(h_ddot_xx(i,k)*z(i,k)*((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+1))/h(i,k)^2+(2*(h_dot_x(i,k))^2*z(i,k)*((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+1))/h(i,k)^3);
        L_zz(i,k)=-lam*((z(i,k)*((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((-(a_23*a_24*h_dot_x(i,k)*h_dot_z(i,k)*(z(i,k)/h(i,k))^a_24)/h(i,k)^2)+(a_23*a_24*h_ddot_xz(i,k)*(z(i,k)/h(i,k))^a_24)/h(i,k)+(a_21*a_22*h_dot_x(i,k)*h_dot_z(i,k)*(z(i,k)/h(i,k))^a_22)/h(i,k)^2-(a_21*a_22*h_ddot_xz(i,k)*(z(i,k)/h(i,k))^a_22)/h(i,k)+(a_23*a_24^2*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_24*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k)-(a_21*a_22^2*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_22*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k))+(2*a_12*h_dot_x(i,k)*h_ddot_xz(i,k)+h_ddot_xz(i,k))*((a_23*a_24*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_24)/h(i,k)-(a_21*a_22*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_22)/h(i,k))+(2*a_12*h_dot_x(i,k)*h_dddot_xxz(i,k)+h_dddot_xxz(i,k)+2*a_12*h_ddot_xz(i,k)*h_ddot_xx(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+(2*a_12*h_dot_x(i,k)*h_ddot_xx(i,k)+h_ddot_xx(i,k))*((a_21*a_22*h(i,k)*(z(i,k)/h(i,k))^a_22*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k)-(a_23*a_24*h(i,k)*(z(i,k)/h(i,k))^a_24*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k))))/h(i,k)-(h_dot_z(i,k)*z(i,k)*((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((a_23*a_24*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_24)/h(i,k)-(a_21*a_22*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_22)/h(i,k))+(2*a_12*h_dot_x(i,k)*h_ddot_xx(i,k)+h_ddot_xx(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)))/h(i,k)^2+((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((a_23*a_24*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_24)/h(i,k)-(a_21*a_22*h_dot_x(i,k)*(z(i,k)/h(i,k))^a_22)/h(i,k))+(2*a_12*h_dot_x(i,k)*h_ddot_xx(i,k)+h_ddot_xx(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20))/h(i,k)-(h_dot_x(i,k)*z(i,k)*((2*a_12*h_dot_x(i,k)*h_ddot_xz(i,k)+h_ddot_xz(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+(a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((a_21*a_22*h(i,k)*(z(i,k)/h(i,k))^a_22*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k)-(a_23*a_24*h(i,k)*(z(i,k)/h(i,k))^a_24*(1/h(i,k)-(h_dot_z(i,k)*z(i,k))/h(i,k)^2))/z(i,k))))/h(i,k)^2+(2*h_dot_x(i,k)*h_dot_z(i,k)*z(i,k)*((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+1))/h(i,k)^3-(h_ddot_xz(i,k)*z(i,k)*((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+1))/h(i,k)^2-(h_dot_x(i,k)*((a_12*(h_dot_x(i,k))^2+h_dot_x(i,k))*((-a_23*(z(i,k)/h(i,k))^a_24)+a_21*(z(i,k)/h(i,k))^a_22+a_20)+1))/h(i,k)^2);
        eps_dot_13(i,k)=(L_xz(i,k)+L_zx(i,k))/2;
        eps_dot_11(i,k)=L_xx(i,k);

        % Calculation of strain values
    if (i==1)

        D_xx(i,k)=0;
        D_xz(i,k)=0;
        D_zx(i,k)=0;
        D_zz(i,k)=0;

        t(i,k)=0;
        dx(i,k)=0;
      else

        D_xx(i,k)=(L_xx(i,k)+L_xx(i-1,k))/(v_x(i,k)+v_x(i-1,k))*(x(i,k)-x(i-1,k))+D_xx(i-1,k);
        D_xz(i,k)=(L_xz(i,k)+L_xz(i-1,k))/(v_x(i,k)+v_x(i-1,k))*(x(i,k)-x(i-1,k))+D_xz(i-1,k);
        D_zx(i,k)=(L_zx(i,k)+L_zx(i-1,k))/(v_x(i,k)+v_x(i-1,k))*(x(i,k)-x(i-1,k))+D_zx(i-1,k);
        D_zz(i,k)=(L_zz(i,k)+L_zz(i-1,k))/(v_x(i,k)+v_x(i-1,k))*(x(i,k)-x(i-1,k))+D_zz(i-1,k);
        D_xx_real(i,k)=abs(L_xx(i,k)+L_xx(i-1,k))/(v_x(i,k)+v_x(i-1,k))*(x(i,k)-x(i-1,k))+D_xx(i-1,k);
        D_xz_real(i,k)=abs(L_xz(i,k)+L_xz(i-1,k))/(v_x(i,k)+v_x(i-1,k))*(x(i,k)-x(i-1,k))+D_xz(i-1,k);
        D_zx_real(i,k)=abs(L_zx(i,k)+L_zx(i-1,k))/(v_x(i,k)+v_x(i-1,k))*(x(i,k)-x(i-1,k))+D_zx(i-1,k);
        D_zz_real(i,k)=abs(L_zz(i,k)+L_zz(i-1,k))/(v_x(i,k)+v_x(i-1,k))*(x(i,k)-x(i-1,k))+D_zz(i-1,k);

        t(i,k)=2/(v_x(i,k)+v_x(i-1,k))*(x(i,k)-x(i-1,k))+t(i-1,k);
        dx(i,k)=lam/h(i,k)*(t(i,1)-t(i,k));
     end

    end
end


% OUTPUT datafile and figures
% save workspace
save FLM_output_data


% Plot the numerical results
fig01=figure (1);
contourf(x,z,z_s);
title("z_s");
xlabel("x");
ylabel("z");
colorbar();

fig02=figure (2);
contourf(x,z,z_n);
title("z_n");
xlabel("x");
ylabel("z");
colorbar();

fig03=figure (3);
contourf(x,z,h);
title("h");
xlabel("x");
ylabel("z");
colorbar();


fig04=figure (4);
contourf(x,z,h_dot_x);
title("h_{dot_x}");
xlabel("x");
ylabel("z");
colorbar();

fig05=figure (5);
contourf(x,z,h_ddot_xx);
title("h_{ddot_xx}");
xlabel("x");
ylabel("z");
colorbar();


fig06=figure (6);
contourf(x,z,v_x);
title("v_x");
xlabel("x");
ylabel("z");
colorbar();

fig07=figure (7);
contourf(x,z,v_z);
title("v_z");
xlabel("x");
ylabel("z");
colorbar();

fig08=figure (8);
contourf(x,z,v_tot);
title("v_{tot}");
xlabel("x");
ylabel("z");
colorbar();

fig09=figure (9);
contourf(x,z,L_xx);
title("L_{xx}");
xlabel("x");
ylabel("z");
colorbar();

fig10=figure (10);
contourf(x,z,L_xz);
title("L_{xz}");
xlabel("x");
ylabel("z");
colorbar();

fig11=figure (11);
contourf(x,z,L_zx);
title("L_{zx}");
xlabel("x");
ylabel("z");
colorbar();

fig12=figure (12);
contourf(x,z,L_zz);
title("L_{zz}");
xlabel("x");
ylabel("z");
colorbar();

fig13=figure (13);
contourf(x,z,D_xx);
title("D_{xx}");
xlabel("x");
ylabel("z");
colorbar();

fig14=figure (14);
contourf(x,z,D_xz);
title("D_{xz}");
xlabel("x");
ylabel("z");
colorbar();

fig15=figure (15);
contourf(x,z,D_zx);
title("D_{zx}");
xlabel("x");
ylabel("z");
colorbar();

fig16=figure (16);
contourf(x,z,D_zz);
title("D_{zz}");
xlabel("x");
ylabel("z");
colorbar();

fig17=figure (17);
contourf(x,z,t);
title("t");
xlabel("x");
ylabel("z");
colorbar();

fig18=figure (18);
contourf(x,z,dx);
title("dx");
xlabel("x");
ylabel("z");
colorbar();

fig19=figure (19);
plot(x,v_x);
title("v_x");
xlabel("x");
ylabel("v_x");

fig20=figure (20);
plot(x,v_z);
title("v_z");
xlabel("x");
ylabel("v_z");

fig21=figure (21);
plot(x,v_tot);
title("v_{tot}");
xlabel("x");
ylabel("v_{tot}");

fig22=figure (22);
plot(x,L_xx);
title("L_{xx}");
xlabel("x");
ylabel("L_{xx}");

fig23=figure (23);
plot(x,L_xz);
title("L_{xz}");
xlabel("x");
ylabel("L_{xz}");

fig24=figure (24);
plot(x,L_zx);
title("L_{zx}");
xlabel("x");
ylabel("L_{zx}");

fig25=figure (25);
plot(x,L_zz);
xlabel("x");
title("L_{zz}");
ylabel("L_{zz}");

fig26=figure (26);
plot(x,eps_dot_13);
xlabel("x");
title("eps_{dot,13}");
ylabel("eps_{dot,13}");

fig27=figure (27);
contourf(x,z,dx);
title("dx");
xlabel("x");
ylabel("z");
colorbar();

fig28=figure (28);
plot(x,dx);
title("dx");
xlabel("x");
ylabel("dx");

fig29=figure (29);
plot(dx(1.25*N,:),z(1.25*N,:));
hold on;
hold off;
title("dx_{final}");
xlabel("dx");
ylabel("z");

fig30=figure (30);
contourf(x,z,eps_dot_11);
title("{\\epsilon}_{11}","interpreter", "tex");
xlabel("x");
ylabel("z");
colorbar();

fig31=figure (31);
contourf(x,z,eps_dot_13);
title("{\\epsilon}_{13}","interpreter", "tex");
xlabel("x");
ylabel("z");
colorbar();


% Save FIGs
print(fig01,'-dpng','FIG01.png');
print(fig02,'-dpng','FIG02.png');
print(fig03,'-dpng','FIG03.png');
print(fig04,'-dpng','FIG04.png');
print(fig05,'-dpng','FIG05.png');
print(fig06,'-dpng','FIG06.png');
print(fig07,'-dpng','FIG07.png');
print(fig08,'-dpng','FIG08.png');
print(fig09,'-dpng','FIG09.png');
print(fig10,'-dpng','FIG10.png');
print(fig11,'-dpng','FIG11.png');
print(fig12,'-dpng','FIG12.png');
print(fig13,'-dpng','FIG13.png');
print(fig14,'-dpng','FIG14.png');
print(fig15,'-dpng','FIG15.png');
print(fig16,'-dpng','FIG16.png');
print(fig17,'-dpng','FIG17.png');
print(fig18,'-dpng','FIG18.png');
print(fig19,'-dpng','FIG19.png');
print(fig20,'-dpng','FIG20.png');
print(fig21,'-dpng','FIG21.png');
print(fig22,'-dpng','FIG22.png');
print(fig23,'-dpng','FIG23.png');
print(fig24,'-dpng','FIG24.png');
print(fig25,'-dpng','FIG25.png');
print(fig26,'-dpng','FIG26.png');
print(fig27,'-dpng','FIG27.png');
print(fig28,'-dpng','FIG28.png');
print(fig29,'-dpng','FIG29.png');
print(fig30,'-dpng','FIG30.png');
print(fig31,'-dpng','FIG31.png');

close all;
