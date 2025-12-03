function [xhistrdf,yhistadjusted]=calculateRDF(particles)
%Calculate the radial distribution function.  The box can be mirrored to
%provide the long range RDF by setting boxx, boxy, boxz to -1:1
rdf=[];
for boxx=0:0
    for boxy=0:0
            %perform every possible subtraction for each of three position
            %variables

            sij(:,:,1)=(bsxfun(@minus,particles.pos(:,1),particles.pos(:,1)')+boxx);
            sij(:,:,2)=(bsxfun(@minus,particles.pos(:,2),particles.pos(:,2)')+boxy);
            %sij(:,:,3)=(bsxfun(@minus,particles.pos(:,3),particles.pos(:,3)')+boxz);
            sij=refold_positions(sij);

            %calculate all of the distances
            rsqij=(sij(:,:,1)*particles.box_size(1)).^2+(sij(:,:,2)*particles.box_size(2)).^2;
            %+(sij(:,:,3)*particles.box_size(3)).^2;
            rij=rsqij.^.5;

            %add the results to the other RDF results (useful for when the
            %box is mirrored)
            rdf=[rdf;rij(:)];
    end
end

%get rid of distance coming from particles and
%themselves
rdf(rdf==0)=[];

%seperate the particle distances into bins
bins=linspace(0,min(particles.box_size/2),100);
bins=0:.02:min(particles.box_size/2)+1;
[yhistrdf xhistrdf]=hist(rdf,bins);

%distances less than half the box size can be trusted,
%ones greater than half cannot
safe=(xhistrdf<(min(particles.box_size)/2)) & xhistrdf~=0;

%Adjust the RDF for the varying volume of the 3D
%spherical shell
for k=2:length(xhistrdf)-1
   yhistadjusted(k)= yhistrdf(k)./4*3/pi./((xhistrdf(k+1)+xhistrdf(k)).^3/8-((xhistrdf(k-1)+xhistrdf(k))).^3/8)./particles.density./particles.n;
end

xhistrdf=xhistrdf(safe);
yhistadjusted=yhistadjusted(safe);
end
