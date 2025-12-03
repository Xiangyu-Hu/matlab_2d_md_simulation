function [xhistrdf yhistave,flag]=md1(in_pos,in_bsize,in_stepnum,in_dt,in_maxt,in_rhorequest,in_temp_request,hand,saveloc)
% General Purpose: This file conducts a simple molecular dynamics
% simulation using lennard jones potentials.  The verlet algorithm is used
% % to solve for the movement of each particle as a result of interparticle
% forces.  As the simulation progresses, the results are shown on a simple
% GUI, as well as well as the progression of several characteristic values.

%Input:
%in_pos(n x 3) - XYZ position of each particle (descending)
%in_bsize - Box size
%in_step_num - total number of steps
%in_dt - time increment for each simulation step
%in_maxt - upper limit on the total length of time
%in_rhorequest - if greater than 0, specify a density (normalizing the box
%size as appropriate).  if less than or equal to zero, calculate the
%density from the given number of particles and the given box size
%in_temp_request - if greater than 0, specify the average particle
%temperature (normalizing the temperatures as necessary).  if less than or
%equal to 0, keep the simulation adiabatic


%Set a necessary variable used later (for the pause button)
cont=false;

%Number of steps between each update of graphs on the GUI
increment=10;

%Get the number of particles from the position vector
num=size(in_pos,1);

%Initialize the necessary properties of the main structure that contains
%information about the simulation
[particles,simcontrol,potential,statistics]=initialize;
simcontrol.step_num=in_stepnum;
simcontrol.deltat=in_dt;

xhist=linspace(0,sum(particles.vel(1,:).^2)^.5*3,40);

veldisty=0;
yhistadjusted=[];
xhistrdf=[];
yhistave=[];
pressuresave=[];
convergence=0.2;
h=0;
fh=hand;
scrsz = get(0,'ScreenSize');
cursize=get(fh,'Position');
cursize(1)=0.4*scrsz(3);cursize(2)=.5*scrsz(4);
set(fh,'Position',cursize);
flag=0;
%call the main function now that all necesary variables are set up
evolve_sample
output_data


    function [particles,simcontrol,potential,statistics]=initialize()

        %Initialize empty structure
        particles=init_particles(num);

        %fill in known fields
        particles.box_size=in_bsize;
        particles.volume=prod(particles.box_size);
        particles.density=particles.n/particles.volume;

        %initialize more empty structures
        simcontrol=init_simulationcontrol;
        statistics=init_statistics;

        %adjust the density if necessary (as described in the input
        %section; in_rhorequest>0)
        simcontrol.rho_requested=in_rhorequest;
        simcontrol.temp_requested=in_temp_request;
        simcontrol.rho_change=simcontrol.rho_requested>0;

        if (simcontrol.rho_change==1)
            scale= (particles.density/simcontrol.rho_requested)^(1/particles.dim);
            particles.box_size=scale*particles.box_size;
            particles.volume=prod(particles.box_size);
            particles.density=particles.n/particles.volume;
        end

        potential=init_potential(particles);

        %assign the particle positions
        particles.pos=in_pos;

        %start with no kinetic/potential energy
        particles.ene_pot(:,:)=0;
        particles.ene_kin(:,:)=0;

        %INSERT CODE FOR GETTING VELOCITIES/ACCEL
        particles.acc(:,:)=0;
        %assign random velocities, and let the calculation routine rescale
        %it later
        particles.vel=2*(rand(size(particles.vel))-.5);
        particles.vel=particles.vel./repmat(sum(particles.vel.^2,2).^.5,[1 particles.dim]);

        temperature=compute_temperature(particles);
        chi=sqrt(simcontrol.temp_requested/temperature);
        particles.vel=chi*particles.vel;

        %calculate the center of mass
        mass_center=sum(particles.pos,1)/particles.n;

        %make sure that the particles' positions are centered around
        %(0,0,0)
        particles.pos=particles.pos-repmat(mass_center,[particles.n 1]);

    end

    function evolve_sample


        %compute temperature
        temperature=compute_temperature(particles);

        %set the sums to 0 (used later for average temperature/pressure)
        temperature_sum=0;
        pressure_sum=0;

        %intialize necessary variables
        pressuresave=zeros(simcontrol.step_num,1);
        cont=true;
        presss=[];
        steps=[];
        virials=[];
        kins=[];
        tots=[];
        yhistadjustedsave=[];
        yhistsavecounter=0;
        velsavecounter=0;
        displsave=0;
        mdisplacement=0;

        %Make sure the positions are within the box
        particles.pos=refold_positions(particles.pos);

        %Build the neighbor list (only done once)
        neighbors=buildneighborslist(particles,potential);

        for step=1:simcontrol.step_num
            %While the 'continue' variable is false, pause the simulation
            %indefinitely
            %             handles()

            fop=fopen('simdata','r');
            stat=fread(fop);
            fclose(fop);
            while stat==1
                fop=fopen('simdata','r');
                stat=fread(fop);
                fclose(fop);
                pause(1)
            end

            if stat==2
                flag=1;
                return;
            end

            %make sure that all particles are in the box
            particles.pos=refold_positions(particles.pos);



            %increment positions based on current accelerations/velocities
            %according to the verlet algorithm
            particles.pos=particles.pos+simcontrol.deltat*particles.vel+.5*simcontrol.deltat^2*particles.acc;

            %increment particle velocities according to the verlet
            %algorithm
            displ=.5*simcontrol.deltat*particles.acc;
            veltemp=particles.vel+.5*simcontrol.deltat*particles.acc;
            
            %Compute forces
            particles=compute_forces(particles,potential,neighbors);
        
            %Rescale velocities for constant temp
            if (simcontrol.temp_constant && temperature>0)
                temperature=compute_temperature(particles);
                chi=sqrt(simcontrol.temp_requested/temperature);
                particles.vel=chi*veltemp+.5*simcontrol.deltat*particles.acc;
            else
                particles.vel=veltemp+.5*simcontrol.deltat*particles.acc;
            end

            %compute temperature, and calculate average potential/kinetic energies
            [temperature,ene_kin_aver,particles]=compute_temperature(particles);
            ene_pot_aver=sum(particles.ene_pot)/particles.n;
            ene_tot_aver=ene_kin_aver+ene_pot_aver;

            %calculate and save the current pressure
            pressure=particles.density*temperature+particles.virial/particles.volume;
            pressuresave(step)=pressure;

            %If it is time to update the graphs, do it
            curpos=ceil(step/increment);
            if mod(step,increment)==1
                presss(curpos)=pressure;
                virials(curpos)=particles.virial;
                kins(curpos)=ene_kin_aver;
                tots(curpos)=ene_tot_aver;
            else
                presss(curpos)=presss(curpos)+pressure;
                virials(curpos)=virials(curpos)+particles.virial;
                kins(curpos)=kins(curpos)+ene_kin_aver;
                tots(curpos)=tots(curpos)+ene_tot_aver;
            end
            [veldistytemp ab]=hist(sum((particles.vel).^2,2).^.5,xhist);
%             veldisty=0;
            veldisty=veldisty+veldistytemp;
            if mod(step,increment)==0

                figure(fh);
                figure(1); 
                %add the current pressure, virial, and step position to the
                %save list
                presss(curpos)=presss(curpos)/increment;
                virials(curpos)=virials(curpos)/increment;
                kins(curpos)=kins(curpos)/increment;
                tots(curpos)=tots(curpos)/increment;
                %                 presss=[presss;pressure];
                steps=[steps;step];
                if step/increment>1
                    cut=ceil((step/increment)*.25);
                    %                                 if curpos>cut
                    p = polyfit(1:cut+1,presss(curpos-cut:curpos),1);
                    err(step/increment)=abs(p(1)/std(presss(curpos-cut:curpos)));
                    %                     disp(sprintf('Ratio of Local Slope to STD: %2.4f',err(step/increment)));
                else
                    err=1; % don't trust results on first pass
                end
                %                                 end
                %                                 presss2=[presss;pressure];
                %                                 virials2=[virials;particles.virial];
                %                                 kins2=[kins;ene_kin_aver];
                %                                 tots2=[tots;ene_tot_aver];

                %plot the position
                                subplot(1,1,1);p1=plot(particles.pos(:,1), particles.pos(:,2) ,'o');xlabel('X');ylabel('Y')
                                xlim([-0.5,0.5])
                                ylim([-0.5,0.5])

                if mod(step,20000)==0
                    1
                end
                yhist=veldisty/sum(veldisty);
                veldisty=0;
                hhold=h;
                htemp=yhist.*log(yhist./(xhist.^2));
%                 htemp=yhist.*log(yhist).*(xhist(2)-xhist(1));
                htemp(isnan(htemp))=[];
                htemp(htemp==Inf)=[];
                htemp(htemp==-Inf)=[];
                h=sum(htemp);
                hhold-h;
                hsave(curpos)=h;

                if  err(step/increment)<convergence
                    thand=title(sprintf('%d steps out of %d (%d%%), Converging ',step,simcontrol.step_num,round(100*step/simcontrol.step_num)));
                else
                    thand=title(sprintf('%d steps out of %d (%d%%), Not Converged',step,simcontrol.step_num,round(100*step/simcontrol.step_num)));
                end

                set(thand,'FontWeight','bold','Fontsize',12);
%                                 set(thand,'Position',[34.36781609195402 64.21487603305785 17.320508075688775],'FontWeight','bold','Fontsize',12);
                if velsavecounter==0 && err(step/increment)<convergence
                    velsave=yhist;
                    velsavecounter=1;
                    yhistsavecounter=0;
                elseif err(step/increment)<convergence
                    velsave=velsave+yhist;
                    velsavecounter=velsavecounter+1;
                    %hold on; plot(xhist,velsave/velsavecounter,'r');hold off
                else
                    velsavecounter=0; %error has gotten worse, reset averages
                end

                %only update the radial distribution function every 4 gui
                %updates
                if mod(length(presss),4)==0
                figure(2);
                %plot the relative distribution of velocities (which should
                %correspond to the boltzmann function at long times)
                vtot=sum((particles.vel).^2,2).^.5;
                bin=0:.02:2;
                [yhist xhist]=hist(sum((particles.vel).^2,2).^.5,20);
                subplot(2,1,1);plot(xhist,yhist,'.');xlabel('Velocity');ylabel('Relative Distribution')
                %subplot(2,3,[1 2]);plot(xhist,yhist,'.');xlabel('Velocity');ylabel('Relative Distribution')
                    %calculate the RDF and trim any interactions coming
                    %from long distances
                    [xhistrdf,yhistadjusted]=calculateRDF(particles);

                    %Plot RDF
                    subplot(2,1,2);plot(xhistrdf,yhistadjusted,'.')
                    if yhistsavecounter==0 && err(step/increment)<convergence
                        yhistadjustedsave=yhistadjusted;
                        yhistsavecounter=yhistsavecounter+1;
                    elseif err(step/increment)<convergence
                        %                                             elseif length(presss)/4>2
                        %Save averages if above 8 times the increment
                        yhistadjustedsave=yhistadjustedsave+yhistadjusted;
                        yhistsavecounter=yhistsavecounter+1;
                        yhistave=yhistadjustedsave/yhistsavecounter;
                        hold on; plot(xhistrdf,yhistadjustedsave/yhistsavecounter,'r');hold off
                    else
                        yhistsavecounter=0; %error has gotten worse, reset averages
                    end
                    xlabel('r')
                    ylabel('g(r)')
                end

                %update the gui
                refresh
                drawnow

                %write current values to the console (just to give an
                %indication of whether the simulation is still stable)
                %disp(sprintf('%3d %0.6f %0.6f %0.6f %0.6f %0.6f',step,temperature,ene_kin_aver,ene_pot_aver,ene_kin_aver+ene_pot_aver,pressure))
            end

            %update the summation variables
            temperature_sum=temperature_sum+temperature;
            statistics.ene_kin_sum=statistics.ene_kin_sum+ene_kin_aver;
            statistics.ene_pot_sum=statistics.ene_pot_sum+ene_pot_aver;
            statistics.pressure_sum=pressure_sum+pressure;

            %Update the neighbor list using the displacements
            [neighbors,mdisplacement]=updateneighborslist(particles,potential,neighbors,displ*max(particles.box_size),mdisplacement);

        end

    end

    function start(hObject,eventdata)
        cont=~cont;
        %         pause(0.5);
    end

    function update_increment(hObject,eventdata)
        %called when the value in the textbox is changed, sets the
        %corresponding variable to the new value
        increment=str2num(get(hObject,'String'));
    end

    function update_dt(hObject,eventdata)
        %called when the value in the textbox is changed, sets the
        %corresponding variable to the new value
        simcontrol.deltat=str2num(get(hObject,'String'));
    end


    function output_data
        dlmwrite(saveloc,'Temperature, Density, Number of particles, Number of time steps, Step size','')
        dlmwrite(saveloc,[in_temp_request in_rhorequest particles.n simcontrol.step_num simcontrol.deltat],'-append');
        dlmwrite(saveloc,'Average potential energy, Pressure','-append','delimiter','')
        dlmwrite(saveloc,[sum(particles.ene_pot)/particles.n mean(pressuresave)],'-append');
        dlmwrite(saveloc,'Radial Distribution Function','-append','delimiter','')
        dlmwrite(saveloc,'r, g(r)','-append','delimiter','')
        dlmwrite(saveloc,[xhistrdf;yhistave]','-append');
    end
end

function [temperature,ene_kin_aver,particles]=compute_temperature(particles)
%compute the temperatures based on the velocities
real_vel=particles.vel.*repmat(particles.box_size,[particles.n 1]);
particles.ene_kin=0.5*sum(real_vel.^2,2);
ene_kin_aver=sum(particles.ene_kin)/particles.n;
temperature=2*ene_kin_aver/particles.dim;
end

function particles=init_particles(n)
%Simple structures corresponding to module: particles
%This structure contains most of the information about the physical
%variables of the simulation
particles.dim=2;
particles.n=n;
particles.box_size=zeros(particles.dim,1);
particles.coord=zeros(particles.n,1);
particles.deltar=zeros(particles.n,1);
particles.density=0;
particles.deru=zeros(particles.n,1);
particles.ene_kin=zeros(particles.n,1);
particles.ene_pot=zeros(particles.n,1);
particles.pos=zeros(particles.n,particles.dim);
particles.vel=zeros(particles.n,particles.dim);
particles.acc=zeros(particles.n,particles.dim);
particles.vel_acc=0;
particles.virial=0;
particles.volume=0;
end

function simcontrol=init_simulationcontrol()
%corresponding to module: simulationcontrol
%This structure contains simulation parameters
simcontrol.deltat=0;
simcontrol.sampin=' ';
simcontrol.sampout=' ';
simcontrol.skin=0;
simcontrol.step_num=0;
simcontrol.temp_constant=1;
simcontrol.temp_requested=100;
simcontrol.title=' ';
simcontrol.rho_requested=0;
simcontrol.rho_change=0;
simcontrol.kb=8.617385D-05;
end

function potential=init_potential(particles)
%corresponding to module: potential
%This structure contains parameters related to the calculation the lennard
%jones interactions
potential.r_cutoff=2.5;
potential.r_min=.5;
potential.buffer=1;
potential.TableSize=20000;
potential.PhiTab=zeros(potential.TableSize,1);
potential.DPhiTab=zeros(potential.TableSize,1);
potential.phicutoff=-100.0/(potential.r_cutoff);
%4/(potential.r_cutoff^12)-4/(potential.r_cutoff^6);
potential.RsqMin=potential.r_min^2;
potential.DeltaRsq=(potential.r_cutoff^2-potential.RsqMin)/(potential.TableSize-1);
potential.InvDeltaRsq=1/potential.DeltaRsq;

%precompute the lennard jones potential in case tabling is needed later
rsqij=.001:.001:max(particles.box_size)^2;
rm2=1./rsqij;
rm6=rm2.*rm2.*rm2; %for some reason, rm2.^3 in matlab is orders of magnitude slower than this
rm12=rm6.^2;
phi=4*(rm12-rm6)-potential.phicutoff;
dphi=24*rm2.*(2*rm12-rm6);
potential.xtable=rsqij;
potential.phitable=phi+dphi*i;
end

function statistics=init_statistics
%corresponding to module: statistics
%This structure contains variables to hold the necessary sums
statistics.ene_kin_sum = 0;
statistics.ene_pot_sum = 0;
statistics.pressure_sum  = 0;
statistics.temperature_sum = 0;
end

function neighbors=buildneighborslist(particles,potential)
%Build a neighbor matrix by examining the distances between each set of
%particles (n^2 operation)

neighbors2=triu(ones(particles.n),1);
neighbors=zeros(particles.n);
[i,j]=find(neighbors2);
inds=sub2ind(size(neighbors),i,j);
diffval=refold_positions(particles.pos(i,:)-particles.pos(j,:));
rsqij=(diffval(:,1)*particles.box_size(1)).^2+(diffval(:,2)*particles.box_size(2)).^2;
%+(diffval(:,3)*particles.box_size(3)).^2;
neighbors(inds)=rsqij<(potential.r_cutoff+potential.buffer)^2;
neighbors=neighbors+neighbors';
end

function [neighbors,mdisplacement]=updateneighborslist(particles,potential,neighbors,displacement,mdisplacement)
%This function updates an existing neighbor list and only considers those
%pairs whose movement could possible change their pair status

%Calculate total distance moved for each particle
displacement=sum(displacement.^2,2).^.5;

%Add each of these distances to the relevant positions in mdisplacement so
%that the maximum possible change for each of the pairs is known
mdisplacement=mdisplacement+bsxfun(@plus,displacement,displacement');

%Find the pairs whose particles have moved enough to warrant an update on
%their neighbor status
combineddist=mdisplacement>potential.buffer;

%The identified pairs will be updated, so their error for the next pass
%will 0
mdisplacement(combineddist)=0;

%Only consider half of the pairs (for symmetry)
neighbors2=triu(combineddist,1);

%Find the nonzero elements, get their indexes (which also represent the
%particles involved), and find the distances for the relevant particles
%pairs
[i,j]=find(neighbors2);
inds=sub2ind(size(neighbors),i,j);
diffval=refold_positions(particles.pos(i,:)-particles.pos(j,:));
rsqij=(diffval(:,1)*particles.box_size(1)).^2+(diffval(:,2)*particles.box_size(2)).^2;
%+(diffval(:,3)*particles.box_size(3)).^2;

%Update the neighbors matrix, and mirror the results (for symmetry)
neighbors(inds)=rsqij<(potential.r_cutoff+potential.buffer)^2;
neighbors(sub2ind(size(neighbors),j,i))=neighbors(inds);

end

