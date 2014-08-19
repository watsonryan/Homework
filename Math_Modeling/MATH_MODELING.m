%%Author: Ryan Watson
%%Date: 8/19/2014
%%Mathematical Modeling Assignment 1
%%
%%
%%Problem 4

A1B1 = [176 , 180];
A1B2 = [176 , 182];
A2B1 = [186 , 182];
A2B2 = [184 , 188];

% Mean
mu_A1B1 = (A1B1(1,1)+A1B1(1,2))/2;
mu_A1B2 = (A1B2(1,1)+A1B2(1,2))/2;
mu_A2B1 = (A2B1(1,1)+A2B1(1,2))/2;
mu_A2B2 = (A2B2(1,1)+A2B2(1,2))/2;

%Standard Dev
std_dev_A1B1 = std(A1B1);
std_dev_A1B2 = std(A1B2);
std_dev_A2B1 = std(A2B1);
std_dev_A2B2 = std(A2B2);

x_A1B1 = -25+mu_A1B1:.01:25+mu_A1B1;
y_A1B1 = (std_dev_A1B1*sqrt(2*pi))^-1*exp(-.5*((x_A1B1-mu_A1B1)/std_dev_A1B1).^2);

x_A1B2 = -25+mu_A1B2:.01:25+mu_A1B2;
y_A1B2 = (std_dev_A1B2*sqrt(2*pi))^-1*exp(-.5*((x_A1B2-mu_A1B2)/std_dev_A1B2).^2);

x_A2B1 = -25+mu_A2B1:.01:25+mu_A2B1;
y_A2B1 = (std_dev_A2B1*sqrt(2*pi))^-1*exp(-.5*((x_A2B1-mu_A2B1)/std_dev_A2B1).^2);

x_A2B2 = -25+mu_A2B2:.01:25+mu_A2B2;
y_A2B2 = (std_dev_A2B2*sqrt(2*pi))^-1*exp(-.5*((x_A2B2-mu_A2B2)/std_dev_A2B2).^2);

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
box(axes1,'on');
hold(axes1,'all');

% Create plot
plot(x_A1B1,y_A1B1,'Parent',axes1,'LineWidth',2,'Color',[1 0 0],...
    'DisplayName','A1 B1');

% Create plot
plot(x_A1B2,y_A1B2,'Parent',axes1,'LineWidth',2,...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'DisplayName','A1 B2');

% Create plot
plot(x_A2B1,y_A1B1,'Parent',axes1,'LineWidth',2,'Color',[0 0 1],...
    'DisplayName','A2 B1');

% Create plot
plot(x_A2B2,y_A1B1,'Parent',axes1,'LineWidth',2,'DisplayName','A2 B2',...
    'Color',[0 0 0]);

% Create xlabel
xlabel({'Mean'});

% Create ylabel
ylabel({'std dev'});

% Create title
title({'Distribution of Shoe/Ball Combinaions'});

% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.799652777777777 0.763726249088668 0.0770833333333333 0.120751988430947]);

% Create textarrow
annotation(figure1,'textarrow',[0.6 0.552604166666667],...
    [0.767672199170125 0.711618257261411],'TextEdgeColor','none',...
    'String',{'mu = 186','sigma = 2.8284'});

% Create textarrow
annotation(figure1,'textarrow',[0.558854166666667 0.5234375],...
    [0.822651452282158 0.780082987551867],'TextEdgeColor','none',...
    'String',{'mu = 184','sigma = 2.8284'});

% Create textarrow
annotation(figure1,'textarrow',[0.367708333333333 0.426041666666667],...
    [0.854809128630705 0.814315352697095],'TextEdgeColor','none',...
    'String',{'mu = 178','sigma = 2.8284'});

% Create textarrow
annotation(figure1,'textarrow',[0.3453125 0.433854166666667],...
    [0.635929460580913 0.569502074688797],'TextEdgeColor','none',...
    'String',{'mu = 179','sigma = 4.2426'});

