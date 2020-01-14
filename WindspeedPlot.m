function [] = WindspeedPlot(W)
figure
max_WS = max(W);
x = 0:max_WS;
y = zeros(length(x),1);
for i = 0:max_WS
    y(i+1) = length(W(W==i)); %calculates the amount of hours for certain windspeeds
end
bar(x,y);                     %creates bar chart with hours (y axis) and windspeed (x axis)
ylabel('Number of hours')     %y axislabel
percent = y/length(W)*100;    %converts the hours into percentage
yyaxis right                  %put y axis to right as percentage

plot(x,percent, 'LineWidth',2)%plots a graph with percent (y axis) and windspeed (x axis)
xlabel('Windspeed')           %x axislabel
ylabel('Percent')             %y axislabel
title('Wind speed distribution')%title
grid on                         %turns on grid
xlim([0 50])                    %x limits
legend('Hours', 'Percent')      %insert legends (first bar chart, 2nd line graph)
end