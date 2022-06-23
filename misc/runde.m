x0 = [1,0.5,2];
t = linspace (0,1000,10000);
y = lsode ("de", x0, t);
# last one is y(end,:)
plot (t,y)