% Create figure of an exponential curve
figure('WindowState','maximized');

X = (1:.001:5);
Y = exp(X);

a = area(Y);
a(1).FaceColor = [0.7 0.7 0.7];
a(1).EdgeColor = [0.7 0.7 0.7];

ax = gca;
box(ax,'off');
set(gca,'visible','off')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'XColor','none')
ax.XAxis.Visible = 'off';