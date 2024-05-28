clearvars
close all

name = 'groundstate'; % file name
load(strcat(name,'.mat'))

%% plot
model = createpde();
Geo = geometryFromMesh(model,T_h.p',T_h.t');

figure(1)

font = 20;
h = gcf;
hold on
clf;
pdeplot(model,"XYData",abs(u_h), "ZData",abs(u_h)','colormap','jet','Mesh','off','FaceAlpha', 1)
view(0,90)
clim([0 1])
ax = gca;
ax.FontSize = font;
map = hot;
map = flipud(map);
colormap(map)
