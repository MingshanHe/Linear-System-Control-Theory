
n = 50;
x = linspace(40,70,n);
y = .7*x + normrnd(0,5,size(x));
plot(x,y)
hold on;
[p,s] = polyfit(x,y,1);
[yfit,dy] = polyconf(p,x,s,'predopt','curve');
fill([x,fliplr(x)],[yfit-dy,fliplr(yfit+dy)],[0.8706 0.9216 0.9804])
line(x,yfit,'color','r')
line(x,yfit-dy,'color','r','linestyle',':')
line(x,yfit+dy,'color','r','linestyle',':')