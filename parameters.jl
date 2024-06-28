#System parameters

N=24 #Horizon

#bounds on inputs
u1low=0;
u1high=100;
u2low=0;
u2high=100;
# u1low=-10;
# u1high=10;
# u2low=-10;
# u2high=10;

#bounds on states
x1low=52.3;
x1high=53.8;
x2low=52.2;
x2high=53.6;
# x1low=-50;
# x1high=53.8;
# x2low=-50;
# x2high=53.6;

# x1low=-100;
# x1high=100;
# x2low=-100;
# x2high=100.0;

# x0=[52.39;52.202];
# t=19;

x0=[53.0,53.0];
xT=[52.7;52.6];
t=0;

A=[0.9867    0.0134; 0.0417    0.9577];
B1=[0.0013    0.0005;0.0008    0.0035];
B2=[-0.0012; -0.0014];
demand=[32.8558, 23.4978, 20.2395, 20.7607, 24.9586, 52.0984, 138.3209, 129.2237, 112.0088, 106.6881, 100.8867, 96.8999, 93.1164, 91.3616, 83.1925, 85.6749, 97.4281, 110.1020, 109.1647, 103.2733, 97.4320, 83.2823, 75.8405, 53.5114]';
# demand=zeros(1,24);
sysC=[0.9956 0;0 0.9793];
sysD=[0.0217 0;0 0.0759];
reservoirPressures=[11;3];
elecPrice=[0.2159, 0.2148, 0.2034, 0.1957, 0.1948, 0.1988, 0.4200, 0.8963, 0.9089, 1.0669, 0.9669, 0.9658, 0.9174, 0.8350, 0.8354, 0.9388, 1.0403, 1.1773, 1.1090, 0.8936, 0.5393, 0.2671, 0.2461, 0.2229]';
sysCb=kron(Matrix{Float64}(I, N, N),sysC);
sysDb=kron(Matrix{Float64}(I, N, N),sysD);
Emodel=[0.054 0 ;0 0.083 ];
Edemand=[0.039 0; 0 0.045];
E=Emodel+Edemand;

# E=[0.0000000001 0;0 0.000000000001];
C=[[1 0;-1 0; 0 1; 0 -1]; zeros(4,2)];
D=[zeros(4,2); [1 0;-1 0; 0 1; 0 -1]];
b=[x1high,-x1low,x2high,-x2low,u1high,-u1low,u2high,-u2low];
Y=[1.0 0;-1 0; 0 1; 0 -1];
z=[x1high, -x1low, x2high ,-x2low];
n=size(A,1);
m=size(B1,2);
l=size(E,2);



# #Big demand setting
# u1low=0;
# u1high=150;
# u2low=0;
# u2high=150;


# #bounds on states
# x1low=52.0;
# x1high=53.8;
# x2low=52.0;
# x2high=53.6;

# Emodel=[0.054 0 ;0 0.083 ];
# Edemand=[0.039 0; 0 0.045];