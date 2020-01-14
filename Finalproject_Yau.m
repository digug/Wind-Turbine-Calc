%ENME 337 Final Project
%Yau, Douglas
%douglas.yau1@ucalgary.ca
clear;clc;close all; %clears variables, command window, and close all matlab windows

%loads all data into matlab for plotting

%loads the windspeed matrix manually created into matlab
windspeed = load('windspeed.mat');
V = struct2cell(windspeed); %converts windspeed from structure into cell
Vdata = cell2mat(V);        %converts from cell to matrix
%take each column of the matrix and stack them into one big column of
%windspeed
V = [Vdata(:,1);Vdata(:,2);Vdata(:,3);Vdata(:,4);Vdata(:,5);Vdata(:,6);Vdata(:,7);Vdata(:,8);Vdata(:,9);Vdata(:,10);Vdata(:,11);Vdata(:,12)];

%loads given chord data that corresonds to radial position
fname = 'chord.dat';        %set the file name as variable fname
fid = fopen(fname, 'r');    %open the data file in read mode
chord = fscanf(fid,'%f');   %scan/read the file into matlab

%loads given twist data that corresonds to radial position
fname = 'twist.dat';        %set the file name as variable fname
fid = fopen(fname, 'r');    %open the data file in read mode
twist = fscanf(fid,'%f');   %scan/read the file into matlab

%loads given radius data for medium blades
fname = 'radius_medium.dat'; %set the file name as variable fname
fid = fopen(fname, 'r');     %open the data file in read mode
radius = fscanf(fid,'%f');   %scan/read the file into matlab

%loads given omega data that corresonds to speeds at hub height
fname = 'omega.dat';        %set the file name as variable fname
fid = fopen(fname, 'r');    %open the data file in read mode
omega = fscanf(fid,'%f');   %scan/read the file into matlab
omega = omega*360/60;

%loads given DU21 data with 3 columns - angle of attack, Cl, Cd respectively
%Cl and Cd corresounds to the angle of attacks
fname = 'DU21.dat';         %set the file name as variable fname
fid = fopen(fname, 'r');    %open the data file in read mode
DU21 = fscanf(fid,'%f');    %scan/read the file into matlab

%define cells for strings of eevry month
months = {'January','February','March','April','May','June','July','August','September','October','November','December'};
%Predifined values provided
rho = 1.23;
ac = 0.2;
Cmin_Wspd = 3;
Cmax_Wspd = 25;
h = 80;
H = 10;
B = 3;
EnergyConsume = 16.5e6;
population = 3519595;   %populations of Montreal in 2016 - from wikipedia
Ptot = 0;               %total power that will be calculated later

AttAngle = DU21(1:3:end);   %find attack angle from DU21
CL = DU21(2:3:end);         %find Cl from DU21
CD = DU21(3:3:end);         %find Cd angle from DU21

V = V(V>Cmin_Wspd & Cmax_Wspd>V);   %get rid of windspeeds that isnt between 3-25


%random values we predefine
a1 = 2;
aprime1 = 2;
comp1 = 1;  %define values > 1e-6 so that we enter the loop at least one
comp2 = 2;
%for the pt matrix - starts at (1,1)
row = 1;
column = 1;

%loops through windspeed 3-25
for i = 3:25
    Vo = i*(h/H)^(1/7);              %Calculates Vo at 80m
    for j = 1:length(radius)                 %loops through the 17 radial positions 
        a = 1;
        aprime = 1;
        while comp1 >1e-6 && comp2 >1e-6     %loops while the difference between a and new a AND aprime and new prime is larger than 1e-6
            psi = atand((1-a)*Vo/((1+aprime)*omega(i)*radius(j)));  %computes psi
            alpha = psi - twist(j);                                 %calculate alpha
            cl = interp1(AttAngle, CL, alpha);                      %interpolate the Cl value at alpha
            cd = interp1(AttAngle, CD, alpha);                      %interpolate the Cd value at alpha
            Cn = cl*cosd(psi)+cd*sind(psi);                         %computes Cn using psi from earlier
            Ct = cl*sind(psi)-cd*cosd(psi);                         %computes Ct using psi from earlier
            sigma = chord(j)*B/(2*pi*radius(j));                    %comptes using radius n chord which is same position as radius
            K = 4*sind((psi))^2/(sigma*Cn);                         %computes K
            a1      = 1/(K+1);                                      %computes new a
            aprime1 = 1/(4*sind(psi)*cosd(psi)/sigma/Ct-1);         %computes new aprime
            if a1>ac                                                %if a is larger than critical value, apply correction
                a1 = .5*(2*K*(1-2*ac)-sqrt((K*(1-2*ac+2))^2+4*(K*ac^2-1)));
            end
            comp1 = abs(a1-a);                                      %compare new a and old a
            comp2 = abs(aprime1-aprime);                            %compare new aprime and old aprime
            a = a1;                                                 %replace old a with new a
            aprime = aprime1;                                       %replace old aprime with new aprime
        end
        comp1 = 1;                                          %reset comparision values so we could enter next loop
        comp2 = 2;                                          %for a calculation and fixing
        Vrel = Vo*(1-a)/sind(psi);                            %Compute Vrel with the fixed a value   
        pt(row,column) = .5*Ct*rho*Vrel^2*chord(j);         %compute pt using Vrel at velocity =i, chord with respect to radius position   

        %the matrix is constructed so that pt value moves down the column
        %with the same speed, but different radial positions.
        %After the column is finished, we move right to the next column
        %size 17x23 - 23 different speeds and 17 radial positions
        %done by using the row and column variable and incrementing them -
        %row resets when column is done (goes back to top)
        row = row+1;                                                
    end
    column = column+1;
    row = 1;
end

%use formula provided to compute integration for M values
for i=1:length(radius)-1
    for j = 1:size(pt,2)
      Ai = (pt(i+1,j)-pt(i,j))/(radius(i+1)-radius(i));
      Bi = (radius(i+1)*pt(i,j)-radius(i)*pt(i+1,j))/(radius(i+1)-radius(i));
      M(i,j) = 1/3*Ai*(radius(i+1)^3-radius(i)^3)+.5*Bi*(radius(i+1)^2-radius(i)^2);
    end
end

Mtot = sum(M)*B;    %sum up each column of M values for total M at each speed
for i = 1:length(Mtot)
    P(i) = Mtot(i).*omega(i+2); %compute power at different speed and omega (omega corresponds to speed)
    if P(i)>5e6
        P(i) = 5e6;
    end
end

%goes through speed data, picks the power that corresponds to the speed and adds to total power
for i = 1:length(V)
    Ptot = Ptot+P(V(i)-2);
end

W = Ptot;
TurbNum = EnergyConsume*population/W;
fprintf('The least number of turbines to satisfy the demands of Montreal is %.0f\n', TurbNum);


for i = 1:12                                %loops for the 12 months of 2016
    figure
    x = Cmin_Wspd:Cmax_Wspd;                %x goes from min cutoff windspeed to max
    y = zeros(length(x),1);                 %Start up y with all zeros (same size as x)
    Vmonth = Vdata(:,i);                    %selects data from corresponding month
    for k = Cmin_Wspd:Cmax_Wspd             %loops through windspeed from 3-25
        y(k-2) = length(Vmonth(Vmonth==k)); %calculates the amount of hours for certain windspeeds
    end
    bar(x,y)                                %plots bar chart for the month
    ylabel('Number of hours')
    
    yyaxis right;
    for j = 1:length(y)
        Pow(j) = P(j)*y(j);                 %plots the power generated at each windspeed
    end
    plot(x,Pow);
    %graph details 
    title(['Windspeed data for ' months{i} ' at 10m'])
    xlabel('Windspeed (m/s)')
    ylabel('Energy (Wh)')
    xlim([Cmin_Wspd-0.6 Cmax_Wspd+0.6]);
    grid on
    legend('Number of hours at speed', 'Energy produced at speed', 'Location', 'eastoutside') 
end