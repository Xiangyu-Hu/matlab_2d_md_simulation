function [particles,potential]=compute_forces(particles,potential,neighbors)

%Reset the quantities that will be calculated
particles.ene_pot(:,:)=0;
particles.acc(:,:)=0;
particles.virial(:,:)=0;

%Get all of the distances between particles.
neighbors2=triu(neighbors,1);
%[i,j,s]=find(neighbors2);
[i,j]=find(neighbors2);
inds=sub2ind(size(neighbors),i,j);
diffval=refold_positions(particles.pos(i,:)-particles.pos(j,:));
rsqij=(diffval(:,1)*particles.box_size(1)).^2+(diffval(:,2)*particles.box_size(2)).^2;
%+(diffval(:,3)*particles.box_size(3)).^2;

%Calculate forces
[phi,dphi]=potential_table(potential,rsqij);

%calculate the potential energy
particles.ene_pot=full(sum(phi));

%calculate the virial
particles.virial=-sum(sum(dphi.*rsqij))/2;

%calcaulate the acceleration for each coordinate (XYZ)
m1=zeros(size(neighbors,1),size(neighbors,2));
m1(inds)=-(dphi.*diffval(:,1));
particles.acc(:,1)=sum(m1,1)-sum(m1,2)';

m1=zeros(size(neighbors,1),size(neighbors,2));
m1(inds)=-(dphi.*diffval(:,2));
particles.acc(:,2)=sum(m1,1)-sum(m1,2)';

%m1=zeros(size(neighbors,1),size(neighbors,2));
%m1(inds)=-(dphi.*diffval(:,3));
%particles.acc(:,3)=sum(m1,1)-sum(m1,2)';

%fix the virial
particles.virial=-particles.virial/particles.dim;

end


function [phi,dphi]=potential_table(potential,rsqij)
rm=1./sqrt(rsqij);
g=1e+2;
rm2=rm.*rm;
phi=-g*rm-potential.phicutoff;
dphi=g*rm2;
end

%%
%function [phi,dphi]=potential_table(potential,rsqij)
%rm2=1./rsqij;
%rm6=rm2.*rm2.*rm2; %for some reason, rm2.^3 in matlab is orders of magnitude slower than this
%rm12=rm6.^2;
%phi=4*(rm12-rm6)-potential.phicutoff;
%dphi=24*rm2.*(2*rm12-rm6);
%end
