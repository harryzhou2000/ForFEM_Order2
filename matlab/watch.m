clear;
load coords
sets = [448                 6832                 6831                  447                 6841                 6839                  455                  457];
plot3(coords(:,2),coords(:,3),coords(:,4),'.',...
    coords(sets,2),coords(sets,3),coords(sets,4),'o')
set(gca,'Clipping','off')