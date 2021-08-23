clear all

T_Final = 1; %Final time in years
T_Step = 4320; %Time step for main simulation in seconds

c = clock;
date = juliandate(c(1), c(2), c(3));
[r1i, v1i] = planetEphemeris(date, 'SolarSystem', 'Earth', '432t', 'km')
[r2i, v2i] = planetEphemeris(date, 'SolarSystem', 'Sun', '432t', 'km')
[r3i, v3i] = planetEphemeris(date, 'SolarSystem', 'Mars', '432t', 'km')
[r4i, v4i] = planetEphemeris(date, 'SolarSystem', 'Venus', '432t', 'km')
[r5i, v5i] = planetEphemeris(date, 'SolarSystem', 'Mercury', '432t', 'km')
[r6i, v6i] = planetEphemeris(date, 'SolarSystem', 'Moon', '432t', 'km')
[r7i, v7i] = planetEphemeris(date, 'SolarSystem', 'Jupiter', '432t', 'km')
[r8i, v8i] = planetEphemeris(date, 'SolarSystem', 'Saturn', '432t', 'km')
[r9i, v9i] = planetEphemeris(date, 'SolarSystem', 'Uranus', '432t', 'km')
[r10i, v10i] = planetEphemeris(date, 'SolarSystem', 'Neptune', '432t', 'km')
[r11i, v11i] = planetEphemeris(date, 'SolarSystem', 'Pluto', '432t', 'km')

names = {'Earth', 'Sun', 'Mars', 'Venus', 'Mercury', 'Moon', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto'}

r1i = r1i';
v1i = v1i';
r2i = r2i';
v2i = v2i';
r3i = r3i';
v3i = v3i';
r4i = r4i';
v4i = v4i';
r5i = r5i';
v5i = v5i';
r6i = r6i';
v6i = v6i';
r7i = r7i';
v7i = v7i';
r8i = r8i';
v8i = v8i';
r9i = r9i';
v9i = v9i';
r10i = r10i';
v10i = v10i';
r11i = r11i';
v11i = v11i';



m1 = 5.972e24;
m2 = 1.989e30;
m3 = .64171e24;
m4 = 4.867e24;
m5 = 3.285e23;
m6 = 7.34767309e22;
m7 = 1.898e27;
m8 = 5.683e26;
m9 = 8.681e25;
m10 = 1.024e26;
m11 = 1.309e22;

stepPlot = false;


m = [5.972e24, 1.989e30, .64171e24, m4, m5, m6, m7, m8, m9, m10, m11]

G = 6.67408e-11;

ui = [r1i; v1i; r2i; v2i; r3i; v3i; r4i; v4i; r5i; v5i; r6i; v6i; r7i; v7i; r8i; v8i; r9i; v9i; r10i; v10i; r11i; v11i].*1000;

G = 6.67408e-11;

T = 43200*365.25*T_Final*2;

dt = T_Step;

disp("Simulating the N-body problem with " + num2str(length(m)) + " bodies for " + num2str(T/86400) + " days with a timestep of " + num2str(dt/86400) + " days")

j = 1;

tolerance = .05;


tic
uf = simulatePlanets(ui, m, dt, T);
toc
plotPlanets(uf, names)

fileName = "PlanetPositions_1y_dt_" + num2str(dt) + ".mat";

save(fileName, 'uf')

[possEclipseInit, possEclipseTime] = findPossEclipses(uf, dt, tolerance);

[eclipseStart, eclipseEnd, type] = simulateEclipses(possEclipseInit, m, possEclipseTime);

data.Eclipse_Start_UTC = eclipseStart';
data.Eclipse_End_UTC = eclipseEnd';
data.Eclipse_Type = type';

dispText = num2str(length(type)) + " Eclipses Found";
disp(dispText);
struct2table(data)

%Funciton that simulates N bodies from initial conditions in ui and masses
%in m, with time step dt and final time T
function [uf] = simulatePlanets(ui, m, dt, T)

u_k = ui;

i = 1;
t = 0;

disp("Simulating Bodies");

while (t<=T)
    
    u_kp1 = u_k + .5*dt*(fcn(u_k, m) + fcn(u_k + dt*(fcn(u_k, m)), m));
    
    uf(:,i) = u_kp1;
    
    u_k = u_kp1;
    
    t = t + dt;
    i = i + 1;
    
end

end


%Function that calculates the gravitational equation for two bodies
function [DVDT] = dvdt(r1, r2, m2, G)
    
    DVDT = (G*m2*(r2-r1))/(norm(r2 - r1))^3;
       
end

%Function to find u dot which is the accelerations and velocities of all
%the bodies being simulated
function [UD] = fcn(u, m)

N = length(m);
j = 1;
UD = [];

G = 6.67408e-11;

DVDTsave = zeros(3,1);

for k = 0:N-1
    DVDT = 0;
    %k;
    for kk = 0:N-1
        if(kk~=k)
            
            DVDT = DVDT + dvdt(u(k*3*2+1:k*3*2+3), u(kk*2*3+1:kk*2*3+3), m(kk+1), G);
            
        end
    end
    
    DVDTsave(:,j) = DVDT;
    j = j + 1;
end

for i = 1:N
    
    UD = [UD; u(i*3*2-2); u(i*3*2-1); u(i*3*2); DVDTsave(:,i)];
    
    
end

end

%Function that creates 2D and 3D plot of simulated system.
function plotPlanets(uf, names)
figure(1)
clf;
hold on;
for ii = 0:length(names)-1
    
    plot(uf(ii*6+1,:), uf(ii*6+2,:), 'linewidth', 2)
    
end

legend(names);

svnm = 'Solar_System1';
print( '-dpng', svnm, '-r200' )

figure(2)
clf;

plot3(uf(1,:), uf(2,:), uf(3,:), 'linewidth', 2)
hold on;

for ii = 1:length(names)-1
    
    plot3(uf(ii*6+1,:), uf(ii*6+2,:), uf(ii*6+3,:), 'linewidth', 2)
    
end

legend(names);

end

function [possEclipseInit, possEclipseTime] = findPossEclipses(uf, dt, tol)
    
    disp("Checking for possible eclipses");
    rE = uf(1:3,:);
    rS = uf(7:9,:);
    rM = uf(31:33,:);
    
    eclipseTime = 0;
    eclipseIndex = 0;
    
    j = 1;
    k = 1;
    
    y = 8115200;
    
    rm = 1737100;
    re = 6378100;
    rs = 696340000;
    
    rSE = rE-rS;
    
    rEM = rM-rE;
    
    possEclipseInit = [];
    possEclipseTime = [];
    

    for i = 1:length(rEM)
        phi = atan((rs + re)/norm(rSE(:,i)));
        
        beta = atan((rs-re)/norm(rSE(:,i)));
        
        x = norm((dot(rSE(:,i), rEM(:,i))/(norm(rSE(:,i))^2))*rSE(:,i));
        
        ymaxe = re - tan(beta)*x + rm/cos(beta);
        
        thetaMaxE = atan(ymaxe/x);
        
        yp = tan(phi)*x + re;
        
        ymaxp = yp + rm/cos(phi);
        
        thetaMaxP = atan(ymaxp/x);
        
        theta = acos((dot(rEM(:,i), rSE(:,i))/(norm(rEM(:,i))*norm(rSE(:,i)))));

        
        if (theta < thetaMaxP + tol && norm(cross(rEM(:,i), rSE(:,i))) > 0 && abs(eclipseTime(k) - i*dt/(86400)) > 5)
            eclipseTime(j) = i*dt/(86400);
            eclipseIndex(j) = i;      
            possEclipseInit(:,j) = uf(:,i-uint32(2*86400/dt));
            possEclipseTime(j) = (i-uint32(2*86400/dt))*dt;
            j = j + 1;
            
            if(j > 2)
                k = k + 1;
            end
        end   

    end
        
                   
    save('eclipseTimes.mat', 'eclipseTime')
    possEclipseInit;
    possEclipseTime;
    eclipseTime;
    eclipseIndex;
    
    dispText = num2str(j-1) + " Possible Eclipses Found";
    disp(dispText)

end

function [startTime, endTime, type] = simulateEclipses(ui, m, possEclipseTime)

T = 86400*4;

dt = 60;

j = 1;

%startTime = [];
%endTime = []
type = "";

for i = 1:length(ui(1,:))
   
    uf = simulatePlanets(ui(:,i), m, dt, T);
    
    dispText = "Simulating Possible Eclipse " + num2str(i) + " of " + num2str(length(ui(1,:)));
    
    disp(dispText);
    
    [startTimeAll(i), endTimeAll(i), typeAll(i)] = findEclipses(uf, dt, possEclipseTime(i));
    
    if(typeAll(i) ~= "No Eclipse")
        startTime(j) = startTimeAll(i);
        endTime(j) = endTimeAll(i);
        type(j) = typeAll(i);
        j = j + 1;
    end
    
end

end

function [startTime, endTime, type] = findEclipses(uf, dt, possEclipseTime)
    
    disp("Checking for eclipses");
    rE = uf(1:3,:);
    rS = uf(7:9,:);
    rM = uf(31:33,:);
    
    y = 8115200;
    
    rm = 1737100;
    re = 6378100;
    rs = 696340000;
    
    rSE = rE-rS;
    
    rEM = rM-rE;
    
    startTime = [ ];
    endTime = [ ];
    type = "";
    
    eclipse = false;
    umbral = false;
    eclipseEnd = false;
    
    startI = 0;
    endI = 0;
    

    for i = 1:length(rEM)
        phi = atan((rs + re)/norm(rSE(:,i)));
        
        beta = atan((rs-re)/norm(rSE(:,i)));
        
        x = norm((dot(rSE(:,i), rEM(:,i))/(norm(rSE(:,i))^2))*rSE(:,i));
        
        ymaxe = re - tan(beta)*x + rm/cos(beta);
        
        thetaMaxE = atan(ymaxe/x);
        
        yp = tan(phi)*x + re;
        
        ymaxp = yp + rm/cos(phi);
        
        thetaMaxP = atan(ymaxp/x);
        
        theta = acos((dot(rEM(:,i), rSE(:,i))/(norm(rEM(:,i))*norm(rSE(:,i)))));

        if(eclipse == false)           
            if (theta < thetaMaxP && norm(cross(rEM(:,i), rSE(:,i))) > 0)
                startI = i;
                eclipse = true;
            end 
        end
        
        if(eclipse == true && eclipseEnd == false)
            if(theta < thetaMaxE && umbral == false)
                umbral = true;
                startI = i;
            end
            if(theta > thetaMaxE && umbral == true)
                endI = i;
                eclipseEnd = true;
            end
            if(theta > thetaMaxP && umbral == false)
                endI = i;
                eclipseEnd = true;
            end
        end
    end
    
    
    if(eclipse == true)
        
        if(umbral == true)
            type = "Umbral";
        else
            type = "Penumbral";
        end
        
        c = clock;
        year = c(1);
        month = c(2);
        day = c(3);
        
        startT = datetime(year, month, day) + seconds(startI*dt) + seconds(possEclipseTime);
        endT = datetime(year, month, day) + seconds(endI*dt) + seconds(possEclipseTime);
    else
        type = "No Eclipse";
        startT = datetime(0, 0, 0);
        endT = datetime(0, 0, 0);
    end
    
    startTime = convertCharsToStrings(datestr(startT));
    endTime = convertCharsToStrings(datestr(endT));
end