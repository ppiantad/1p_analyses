x7 = s7pcaData(1,:) 
y7 = s7pcaData(2,:) 
z7 = s7pcaData(3,:) 

smoothX7 = imgaussfilt(x7,2) 
smoothY7 = imgaussfilt(y7,2)
smoothZ7 = imgaussfilt(z7,2)


%we can split this into the "avoid period (-.5 to 1s around avoid and the
%rest of the trial around that) 

preAvX7 = smoothX7(1,1:16)
preAvY7 = smoothY7(1,1:16)
preAvZ7 = smoothZ7(1,1:16)
avX7 = smoothX7(1,16:31)
avY7 = smoothY7(1,16:31) 
avZ7 = smoothZ7(1,16:31) 
postAvX7 = smoothX7(1,31:41); 
postAvY7 = smoothY7(1,31:41); 
postAvZ7 = smoothZ7(1,31:41); 
%plot3(x7,y7,z7) 
plot3(preAvX7,preAvY7,preAvZ7,'Color','b','LineStyle','-','LineWidth',2); 
hold on 
plot3(avX7,avY7,avZ7,'Color','k','LineWidth',2)
plot3(postAvX7,postAvY7,postAvZ7,'Color','b','LineStyle','-','LineWidth',2);


%plot3(smoothX7,smoothY7,smoothZ7)
grid on 
hold on


xlabel('PC1', 'FontName', 'Franklin Gothic Medium');
ylabel('PC2' ,'FontName','Franklin Gothic Medium');
zlabel('PC3','FontName','Franklin Gothic Medium')

ax.XAxis.FontName = 'Franklin Gothic Medium';
ax.YAxis.FontName = 'Franklin Gothic Medium';
ax.ZAxis.FontName = 'Franklin Gothic Medium';
%% 
x1 = projS1(1,:) 
y1 = projS1(2,:) 
z1 = projS1(3,:) 

smoothX1 = imgaussfilt(x1,2)
smoothY1 = imgaussfilt(y1,2) 
smoothZ1 = imgaussfilt(z1,2) 
%plot3(x1,y1,z1)
preAvX1 = smoothX1(1,1:16)
preAvY1 = smoothY1(1,1:16)
preAvZ1 = smoothZ1(1,1:16)
avX1 = smoothX1(1,16:31)
avY1 = smoothY1(1,16:31) 
avZ1 = smoothZ1(1,16:31) 
postAvX1 = smoothX1(1,31:41); 
postAvY1 = smoothY1(1,31:41); 
postAvZ1 = smoothZ1(1,31:41); 
%plot3(x7,y7,z7) 
plot3(preAvX1,preAvY1,preAvZ1,'Color','b','LineStyle','-','LineWidth',2); 
hold on 
plot3(avX1,avY1,avZ1,'Color','g','LineWidth',2)
plot3(postAvX1,postAvY1,postAvZ1,'Color','b','LineStyle','-','LineWidth',2);

%% 
x2 = projS2(1,:) 
y2 = projS2(2,:) 
z2 = projS2(3,:) 

smoothX2 = imgaussfilt(x2,2)
smoothY2 = imgaussfilt(y2,2) 
smoothZ2 = imgaussfilt(z2,2) 

plot3(smoothX2,smoothY2,smoothZ2)

%% 
x3 = projS3(1,:) 
y3 = projS3(2,:) 
z3 = projS3(3,:) 

smoothX3 = imgaussfilt(x3,2)
smoothY3 = imgaussfilt(y3,2) 
smoothZ3 = imgaussfilt(z3,2) 

plot3(smoothX3,smoothY3,smoothZ3)

%% 
x4 = projS4(1,:) 
y4 = projS4(2,:) 
z4 = projS4(3,:) 

smoothX4= imgaussfilt(x4,2)
smoothY4 = imgaussfilt(y4,2) 
smoothZ4 = imgaussfilt(z4,2) 
%plot3(x1,y1,z1)
plot3(smoothX4,smoothY4,smoothZ4)

%% 
x5 = projS5(1,:) 
y5 = projS5(2,:) 
z5 = projS5(3,:) 

smoothX5 = imgaussfilt(x5,2)
smoothY5 = imgaussfilt(y5,2) 
smoothZ5 = imgaussfilt(z5,2) 
%plot3(x1,y1,z1)
plot3(smoothX5,smoothY5,smoothZ5)

%% 
x6 = projS6(1,:) 
y6 = projS6(2,:) 
z6 = projS6(3,:) 

smoothX6 = imgaussfilt(x6,2)
smoothY6 = imgaussfilt(y6,2) 
smoothZ6 = imgaussfilt(z6,2) 
%plot3(x1,y1,z1)
preAvX6 = smoothX6(1,1:16)
preAvY6 = smoothY6(1,1:16)
preAvZ6 = smoothZ6(1,1:16)
avX6 = smoothX6(1,16:31)
avY6 = smoothY6(1,16:31) 
avZ6 = smoothZ6(1,16:31) 
postAvX6 = smoothX6(1,31:41); 
postAvY6 = smoothY6(1,31:41); 
postAvZ6 = smoothZ6(1,31:41); 
%plot3(x7,y7,z7) 
plot3(preAvX6,preAvY6,preAvZ6,'Color','b','LineStyle','-','LineWidth',2); 
hold on 
plot3(avX6,avY6,avZ6,'Color','r','LineWidth',2)
plot3(postAvX6,postAvY6,postAvZ6,'Color','b','LineStyle','-','LineWidth',2);


%% 
scatter3(smoothX1(1,21),smoothY1(1,21),smoothZ1(1,21))
scatter3(smoothX7(1,21),smoothY7(1,21),smoothZ7(1,21))
scatter3(smoothX6(1,21),smoothY6(1,21),smoothZ6(1,21))

