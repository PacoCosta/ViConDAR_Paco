as=find(ismembertol(Output.TS.fullWF.time,Output.TS.lidar.time{1,1}));
mean_error_vel_probe=mean(abs(Output.TS.fullWF.Uval{1,1}(1,as)-Output.TS.lidar.Uval{1,1}))
