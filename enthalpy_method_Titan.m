

%% setting up the domain
nx=200; %x nodes
ny=200; %y nodes
p=ceil(ny/2); %middle of y nodes
dx=.105;
dy=dx;
%% Rab conversion
l=10000; %length scale (10000 m)
Rabi= Rab*(width/(dx*l)); 
%% Input parameter values (Titan)
beta = 0.0013; %basement slope
q0= Rabi*beta; %initial sediment flux (dimensionless q in eq.13)
w = 0; %width of sediment strip/2 - 1
nshore = 60; %initial shoreline
B=2*pi/period; 
round_c = 10; %rounding of floating-point numbers（4=Earth, 10=Titan）
%% Time vector
dt=.005;
t=0:dt:Tmax*period;
nt=length(t);
%% Initialize time vectors
s=zeros(1,nt); % Shoreline position
r=zeros(1,nt); % Shoreline position
V = zeros(1,nt);
Zp = zeros(1,nt);
%% Initialize space vectors
%cycle values
x=zeros(1,nx);
y=zeros(1,ny);

FN=zeros(nx,ny); % Sediment flux on northern boundary
FS=zeros(nx,ny); % Sediment flux on southern boundary
FE=zeros(nx,ny); % Sediment flux on eastern boundary
FW=zeros(nx,ny); % Sediment flux western boundary
h=zeros(nx,ny); % River bed elevation respect to SL
H=zeros(nx,ny); % Sediment thickness 
Hnew=zeros(nx,ny); % Sediment thickness 
E=zeros(nx,ny); % Bedrock Elevation 

for i=1:nx
    x(i)=(i-nshore)*dx;
end

for j=1:ny
    y(j)=(j-p)*dy;
end

for i=1:nshore
    E(i,:) = -beta*x(i);
    h(i,:) = -beta*x(i);
end

for i=nshore+1:nx
    E(i,:) = -beta*x(i);
end


FN(1,p-w:p+w) = q0; %Strip of sediment flux

Z = 0;
for k = 1:nt
    
    %%Sea Level Rise
    Z0 = Z;
    Z = A*sin(B*(k*dt));

    dz = Z - Z0;
    Zp(k) = Z;
    
    %% Rest of land domain
    for j = 1:ny
        for i = 1:nx-1
            
            if j ~= 1
                FW(i,j) = (h(i,j) - h(i,j-1))/dx;
            end
            
            if j ~= ny
                FE(i,j) = (h(i,j+1) - h(i,j))/dx;
            end
            
            FS(i,j) = (h(i,j) - h(i+1,j))/dx;
            
            FS(i,j) = min(H(i,j)*dx/dt + FN(i,j) + FE(i,j) - FW(i,j),FS(i,j));

            if FS(i,j) < 0.
                FS(i,j) = 0.;
            end

            FN(i+1,j) = FS(i,j);

        end
    end

    Hnew(:,:) = H(:,:) + dt/dx*(FN(:,:) + FE(:,:) - FS(:,:) - FW(:,:));
    Hnew(:,:) = round(Hnew,round_c); %rounding of floating-point numbers
    H(:,:) = Hnew(:,:);
    h(:,:) = max(H(:,:)+E(:,:)-Z,0);

    if mod(k,1/dt) == 0
        disp(['t=', num2str(k*dt)]);
    end

    %% Seed new values
    H = Hnew;

    if mod(k*dt, period*0.25) == 0
       outfile = ['output/', num2str(scenario),'_tstep_', num2str(dt) ,'_' , num2str(k*dt/period) , '.mat'];
       save(outfile)
    end
end

