[X,Y,Z]=meshgrid(windfield.grid.t, gridy,gridz);
figure,
surf(Y,X,Z,contour (windfield.grid.y,windfield.grid.z,(squeeze(windfield.u(:,1,:)))','Fill','on'))



[z y] = meshgrid(-90:4.5:90);
figure,

hold on,
for i=230.625:10*1.875:270

x = zeros(41,41)+i;



surf(x,y,z);
end
set(gca,'xlim', [225,275])