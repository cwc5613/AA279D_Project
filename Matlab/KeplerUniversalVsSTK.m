function KeplerUniversalVsSTK()
%% Constants
mu = 398600.4418;

r0 = repmat(PosVel(1,1:3)',[1 size(PosVel,1)]);
v0 = repmat(PosVel(1,4:6)',[1 size(PosVel,1)]);
[r,v] = keplerUniversal(r0,v0,0:60:(size(PosVel,1)-1)*60,mu);
PlotStateVectorErrors(r,v,PosVel,0:60:(size(PosVel,1)-1)*60,'Parabolic Orbit Propagation');

end

function PlotStateVectorErrors(r,v,PosVel,t,titleText)
MaxPositionError = abs(r-PosVel(:,1:3)');
MaxVelocityError = abs(v-PosVel(:,4:6)');
figure('color',[1 1 1]);
subplot(4,1,1); plot(t./60,MaxPositionError); title(titleText); ylabel('Pos Err (km)');
subplot(4,1,2); plot(t./60,MaxVelocityError); ylabel('Vel Err (km/s)'); 
subplot(4,1,3); plot(t./60,abs(sqrt(sum(r.^2,1)) - sqrt(sum(PosVel(:,1:3)'.^2,1)))); ylabel('Pos Mag Err (km)'); xlabel('Time (m)');
subplot(4,1,4); plot(t./60,abs(sqrt(sum(v.^2,1)) - sqrt(sum(PosVel(:,4:6)'.^2,1)))); ylabel('Vel Mag Err (km/s)'); xlabel('Time (m)');
end