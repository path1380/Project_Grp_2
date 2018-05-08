clear
nx = 201;
ny = 201;

vidObj = VideoWriter('delta.avi');
open(vidObj);

a1 = dir('sol01_*.txt');
u = load(a1(1).name);
r0 = reshape(u(:,3),nx,ny);
x = reshape(u(:,1),nx,ny);
y = reshape(u(:,2),nx,ny);
r = load(a1(1).name);
r = reshape(r(:,3),nx,ny);

% surf(x,y,r);
% axis tight
% set(gca,'nextplot','replacechildren');

for i = 1:3:length(a1)
    r = load(a1(i).name);
    r = reshape(r(:,3),nx,ny);
%     contour(x,y,r,linspace(-0.5,0.5,40))
    contour(x,y,r,50)
%     colorbar
%      surf(x,y,r)
     colorbar
    %view(90,0)
    currFrame = getframe(gcf,[74 47 450 350]);
    writeVideo(vidObj,currFrame);
end
close(vidObj);