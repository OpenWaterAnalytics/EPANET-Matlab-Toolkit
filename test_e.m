inpname='bwcn8-168-true-ec.inp';
d=epanet(inpname);
% delete controls and add them back
controls = d.getControls;
d.deleteControls;
d.addControls(controls)
controls1 = d.getControls;