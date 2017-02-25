% Improved CPU-based metropolis algorithm for simulating the Ising spin
% lattice; ramps temperature up from 10% below the critical temp to 10%
% above it.

doprint = 1;
%domovie = 1;
Tc = 2.27;              % Critical temperature @ 2.27 K
Ti = Tc;                % Initial temperature for simulation
Tf = Tc;                 % Final temperature for simulation
size = 400;              % Lattice size
numGenerations = 31250;  % Number of times to evaluate the entire lattice.
                        % It's possible not every site will be evaluated.
logicalLattice = (rand(size,size) > .5);  % Random initial lattice, 0s/1s.
initialLattice = 2.*(logicalLattice+1)-3; % Initial lattice, -1s/1s.
magnetization = sum(sum(initialLattice)); % This will represent the 
                                          % magnetization of a given
                                          % generation.
%% Initial plotting of stuff
subplot(1,2,1); % Lattice display on the left
colormap([1 1 1;0 0 0]); % Colors (black and white)
% Information on %done, etc
plot_title = ['Generation 0 out of ' num2str(numGenerations)];
title(plot_title);
imagesc(initialLattice); % Displays the initial lattice
% Saves the current generation frame:
subplot(1,2,2); % Magnetization display on the right
plot(1:1,magnetization);
title('Overall magnetization');
if doprint == 1, print(strcat('generation_',num2str(0)),'-dpng','-r0'); end
%if domovie == 1, mov(1) = getframe; end
drawnow;
genTimer = 0;

%% Evaluate the entire lattice for all of the generations
lattice = initialLattice;
%avgGenTimer = 0;
%avgSiteTimer = 0;
T_step = (Tf-Ti)/numGenerations;
T = Ti;

for gen = 1:numGenerations % For every generation:
    timer = tic(); % See how long each generation takes.
    magnetization(1,gen+1) = sum(sum(lattice)); % Sum the overall magnetization
    T = 2.27 + 2.26*sin(gen*.001);
    
    % Lattice printing stuff
    printLattice = lattice > 0; % Turns the -1 and 1 lattice back to 0s
                                % and 1s
    subplot(1,2,1);
    colormap([1 1 1; 0 0 0]);
    format short
    plot_title = ['Generation ' num2str(gen) ...
        ' out of ' num2str(numGenerations) ' at T=' ...
        num2str(T) '.'];
    imagesc(printLattice);
    title(plot_title);
    xlabelstuff = ['avg time per generation: ' num2str(genTimer) ...
        ' seconds'];
    xlabel(xlabelstuff);
    subplot(1,2,2);
    plot(1:gen+1,magnetization);
    %axis([0 numGenerations -size^2 size^2]);
    title('Overall magnetization');
    drawnow;
    if doprint == 1, print(strcat('generation_',num2str(gen)),'-dpng','-r0'); end
    %if domovie == 1, mov(gen+1) = getframe; end
    
    % Primary analysis for loop per generation
    for latIter = 1:size^2 % For every site:
%         timer = tic();
        x = randi(size,1,1); % Pick a random column.
        y = randi(size,1,1); % pick a random row.
        energyDiff = deltaU(x,y,lattice,size); % Calculate the energy
                                               % should the site be
                                               % flipped.
        
        % Flip lattice site if energy is less than 0, or probabilistically
        % based on the partition function.
        if energyDiff <= 0
            lattice(y,x) = -lattice(y,x);
        else
            if rand < exp(-2*energyDiff/T)
                lattice(y,x) = -lattice(y,x);
            end
        end
%         siteTimer = toc(timer);
%         avgSiteTimer = (siteTimer + avgSiteTimer)/latIter;
    end
    genTimer = toc(timer);
end
% if domovie == 1, VideoWriter(mov,'Ising_lattice_movie','Compression','Cinepak'); end
%% Primary energy state analysis function for energy difference of a
%  potential spin flip on the lattice.
function dU = deltaU(x,y,lattice,size)
    % If the point is at the top of the lattice, then the site above it
    % is at the bottom of the lattice (the world is round).
    if y == 1, top = lattice(size,x); else, top = lattice(y-1,x); end
    % Same, for a point at the bottom.
    if y == size, bottom = lattice(1,x); else, bottom = lattice(y+1,x); end
    % For the left.
    if x == 1, left = lattice(y,size); else, left = lattice(y,x-1); end
    % For the right.
    if x == size, right = lattice(y,1); else, right = lattice(y,x+1); end
    
    % Calculate energy difference for that site.
    dU = 2*lattice(y,x)*(top + bottom + right + left);
end