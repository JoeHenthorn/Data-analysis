clear all, close all

figure(1) = figure('Position', [0.1, 0.1, 1400, 400]);
subplot(1,5,1)
x = linspace(-10,10,1000);
y = 1./(1+exp(-x));
plot([0 0],[0 1],'k'), hold on
plot(x,y,'k','Linewidth',3)
grid on
title({'Sigmoid function'})
xlabel('Input');
ylabel('Activation');
set(gca,'Fontsize',[12])


subplot(1,5,2)
x = linspace(-10,10,1000);
y = tanh(x);
plot([-10 10],[0 0],'k'), hold on
plot([0 0],[-1 1],'k'), hold on
plot(x,y,'k','Linewidth',3)
grid on
title({'Hyperbolic tangent (tanh) function'})
xlabel('Input');
ylabel('Activation');
set(gca,'Fontsize',[12])

% ReLU function
subplot(1,5,3)
x = linspace(-10,10,1000);
y = poslin(x); % Positive Linear transfer function
plot([-10 10],[0 0],'k'), hold on
plot([0 0],[0 10],'k'), hold on
plot(x,y,'k','Linewidth',3)
grid on
title({'ReLU function'})
xlabel('Input');
ylabel('Activation');
set(gca,'Fontsize',[12])

subplot(1,5,4)
x1 = linspace(-10,0,500);
x2 = linspace(0,10,500);
y1 = 0.1*x1; % ?*x
y2 = poslin(x2); % Positive Linear transfer function
y = [y1 y2];
plot([-10 10],[0 0],'k'), hold on
plot([0 0],[-2 10],'k'), hold on
plot(x,y,'k','Linewidth',3)
grid on
title({'Leaky ReLU function'})
xlabel('Input');
ylabel('Activation');
set(gca,'Fontsize',[12])


subplot(1,5,5)
x1 = linspace(-10,0,500);
x2 = linspace(0,10,500);
y1 = x2;
y2 = (exp(x1)-1); % ?*(e^x -1)
y = [y2 y1];
plot([-10 10],[0 0],'k'), hold on
plot([0 0],[-2 10],'k'), hold on
plot(x,y,'k','Linewidth',3)
grid on
title({'ELU function'})
xlabel('Input');
ylabel('Activation');
set(gca,'Fontsize',[12])

saveas(gcf,'Activation function plots.png')
