
[velocity, ~, ~, ~, ~, dz, dy, dt, ~, ~, SummVars]=readBLgrid('C:\Users\fcosta\SWE_LOCAL\ViConDAR_Paco\OriginalWF\QlunDAR_Sh00_SD10_V15_TI25.wnd');
windfield=velocity2windfield(velocity,dz,dy,dt,SummVars);
save('C:\Users\fcosta\SWE_LOCAL\ViConDAR_Paco\OriginalWF\QlunDAR_TI25\QlunDAR_Sh00_SD10_V15_TI25.mat','windfield');
