clear all
clc
%This example is a simple 2D aerofoil in unsteady flow modelled using lumped vortex elements. 
%The singularities associated with ideal irrotational vortices are bypassed by the introduction of the Lamb-Oseen vortex core model. 

%% Define Variables
Tstep = input('Number of time steps = ');  %Number of time steps
vortic = zeros(Tstep,3);         %Declare array to store wake data
rho = 1;                         %Air Density
ut = input('Vx = ');             %Velocity in X-component of the aerofoil
c = 1.56/ut;                     %Chord
h0 = c*0.019;                    %Heaving amplitude
omega = 2*ut/c*8.5;     %Heaving frequency || %If omega = 0, this is just the aerofoil
beta = 0*5*pi/180;    %Orientation of blade with respect to Global Coordinate System
sn = sin(beta);                  %Sine Alpha
cs = cos(beta);                  %Cosine Alpha
dt = 0.0009*c/ut;                %Time step
t = 0:dt:(Tstep-1)*dt;           %Time vector
wt = h0*omega*cos(omega*t); %Z-component of the velocity of the aerofoil

%% Define Location 
%Origin Location of the Body Coordinate System
sx = -ut*t;
sz = -h0*sin(omega*t);

%Location of most recently shed wake element
Q = sqrt(ut^2 + wt.^2);
dxw = 0.3*Q*dt;

%Trailing vortex locations
vortic(:,1) = (c + dxw)*cs + sx;
vortic(:,2) = -(c + dxw)*sn + sz;

%Influence matrix terms
a = -1/(pi*c);
b = 1./(2*pi*(c/4 + dxw));
xx1 = (0.75*c*cs + sx(2:Tstep))';
zz1 = (0.75*c*sn + sz(2:Tstep))';

%Location of the bound vortex
xx1_1 = 0.25*c*cs + sx;
zz1_1 = 0.25*c*sn + sz;

%RHS1_a:
RHS1_a = -(ut*sn + wt*sn);

%Vortex core for simple representation
cutoff = 0.001;

%% Define Arrays:
t1 = zeros(Tstep,1);
L = t1;
D = t1;    
cl = t1;
cd = t1;
clt = t1;
sx1 = t1;

%% Looping wrt Time steps
for it = 1:Tstep
    rhs2 = 0;
    wwake = 0;
    if it == 1      
    else
        %Calculate the induced velocity arising from all the wake elements.
        it1 = it - 1;   %Index marker
        
        AoA_x_pos = 0.75*c*cs + sx(it); %3/4 chord point
        AoA_z_pos = 0.75*c*sn + sz(it); %3/4 chord point
        
        %Distance from wake vortex elements to 3/4 chord point
        rx = AoA_x_pos - vortic(1:it1,1);
        rz = AoA_z_pos - vortic(1:it1,2);
        r2 = rx.^2 + rz.^2 + 0*cutoff^2;
        
        %Induced velocity by each vortex wake element
        u1 = rz/(2*pi).*vortic(1:it1,3)./r2.*(1 - exp(-r2/cutoff^2));
        w1 = -rx/(2*pi).*vortic(1:it1,3)./r2.*(1 - exp(-r2/cutoff^2));        

        %Get rid of NaN values due to singularity at core
        u1(isnan(u1)) = 0;
        w1(isnan(w1)) = 0;   
        
%This needs to be done because in spite of the vortex core model, the code will still get confused at the core's centre. 
%The programme is trying to evaluate 0*0/0, which returns NaN. So we apply to remove it. 
        
        %Sum all velocity contributions from the wake vortex elements
        u = sum(u1);
        w = sum(w1);
        
        %Dot induced velocity with normal vector of aerofoil
        wwake = u*sn + w*cs;
        rhs2 = -sum(vortic(1:it1,3));    
    end
    
    %Equations 
    rhs1 = -ut*sn - wwake - wt(it)*cs;
    vortic(it,3) = 1/(b(it)/a - 1)*(rhs1/a-rhs2);
    gammat = rhs2 - vortic(it,3); 
    
    if it < 1  
    else        
%Velocity induced by the bound vortex on the aerofoil at the positions of all the wake elements
        Bound_x = 0.25*c*cs + sx(it);
        Bound_z = 0.25*c*sn + sz(it);        
        
%Distance from wake vortex elements to aerofoil bound vortex
        rx = vortic(1:it,1) - Bound_x;
        rz = vortic(1:it,2) - Bound_z;        
        r2 = rx.^2 + rz.^2 + 0*cutoff^2;
        
%Induced velocity on each vortex wake element by bound vortex
        u = rz/(2*pi)*gammat./r2.*(1 - exp(-r2/cutoff^2));
        w = -rx/(2*pi)*gammat./r2.*(1 - exp(-r2/cutoff^2));

%Get rid of NaN values due to singularity at core
        u(isnan(u)) = 0;
        w(isnan(w)) = 0;        
        
%Determine the influence of the wake elements on each other
        for i = 1:it

%Distance between all wake elements and the wake element of interest
            rx = vortic(i,1) - vortic(1:it,1);
            rz = vortic(i,2) - vortic(1:it,2);
            r2 = rx.^2 + rz.^2 + 0*cutoff^2;
            
            u11 = rz/(2*pi).*vortic(1:it,3)./r2.*(1 - exp(-r2/cutoff^2));
            w11 = -rx/(2*pi).*vortic(1:it,3)./r2.*(1 - exp(-r2/cutoff^2));
            
%Get rid of NaN values due to singularity at core
            u11(isnan(u11)) = 0;
            w11(isnan(w11)) = 0;
            
%Sum influence of all wake elements on the element of interest
            u1 = sum(u11);
            w1 = sum(w11);
            
            u(i) = u(i) + u1;
            w(i) = w(i) + w1;
            %uw(i,1) = vortic(i,1) + u(i)*dt;
            %uw(i,2) = vortic(i,2) + w(i)*dt;      
        end
        
%Move the wake elements with the local velocities
        vortic(1:it,1) = vortic(1:it,1) + u(1:it)*dt;
        vortic(1:it,2) = vortic(1:it,2) + w(1:it)*dt;
    end
    
    if it == 1
        %No circulation at t = 0
        gamat1 = 0;
    end
    
    %Kinematic Velocity
    kin_vel = 0.5*rho*(ut^2 + wt(it)^2);

    %Rate of change of bound circulation with respect to time
    dgamdt = (gammat - gamat1)/dt;
    
    %Store current value of bound circulation for use in the unsteady Bernoulli equation at the next time step
    gamat1 = gammat;
    AoA_x_pos = 0.75*c*cs + sx(it);
    AoA_z_pos = 0.75*c*sn + sz(it);
    
    %Calculate Downwash
    rx = AoA_x_pos - vortic(1:it,1);
    rz = AoA_z_pos - vortic(1:it,2);
    r2 = rx.^2 + rz.^2 + cutoff^2;
    u11 = rz/(2*pi).*vortic(1:it,3)./r2;
    w11 = -rx/(2*pi).*vortic(1:it,3)./r2;
            
    %Sum influence of all wake elements on the element of interest
    u = sum(u11);
    w = sum(w11);
    
    %Wake-induced Downwash
    ww = u*sn + w*cs;
    
    %Lift Force
    L(it) = rho*(ut*gammat + dgamdt*c);
    
    %Drag Force
    D(it) = rho*(-ww*gammat + dgamdt*c*sn);
    
    %Lift Coefficient
    cl(it) = L(it)/kin_vel/c;
    
    %Drag Coefficient
    cd(it) = D(it)/kin_vel/c;
    
%Note the drag force is a lift-induced drag. It is because of the vorticity of the wake.
    
    sx1(it) = sx(it) - ut*dt;
end

%% Plotting
%First_Plot
if dt>0.001*c/ut
    figure
    plot(vortic(:,1)/c,vortic(:,2)/c,'--')
else
    figure
    plot(vortic(:,1)/c,vortic(:,2)/c,'x')
end

%Plot_1 Labels
xlabel('x/c')
ylabel('z/c')
title('Xcomponent Chord wrt Z-component Chord')
grid on;

%Second_Plot
if omega == 0
    ylim([-0.75 0.75]) %Set Limit
end
figure;
plot(ut*t/c,cl)
grid on;

%Plot_2 Labels
xlabel('u_\infty t/c')
ylabel('C_l')
title('Lift Coefficient wrt X-component Velocity of the aerofoil')
figure;
plot(ut*t/c,cd)

%Plot_3 Labels
xlabel('u_\infty t/c')
ylabel('C_d')
title('Drag Coefficient wrt X-component Velocity of the aerofoil')
grid on; axis equal;
