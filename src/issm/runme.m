steps=[1:3];
set_paths;
addpath('../data/')

if any(steps==1) 
    disp('	Step 1: Mesh');

    % Load a GlaDS-MATLAB mesh
    dmesh = load('../data/mesh.mat');

    md = model;
    md.mesh = mesh2d();
    md.mesh.x = dmesh.tri.nodes(:, 1);
    md.mesh.y = dmesh.tri.nodes(:, 2);
    md.mesh.elements = dmesh.tri.connect;
    md = meshconvert(md,md.mesh.elements,md.mesh.x,md.mesh.y);

    save('MassCon.mesh.mat', 'md')
end 

if any(steps==2) 
    disp('	Step 2: Parameterization');
    md=loadmodel('MassCon.mesh.mat');

    md=setmask(md,'','');

    % Run parameterization script to set up geometry, velocity, material properties, etc.
    md=parameterize(md,'MassCon.par');

    % GLADS HYDROLOGY PARAMETERIZATION
    md.hydrology=hydrologyglads();

    % PARAMETERS));
    md.hydrology.sheet_conductivity = 0.05*ones(md.mesh.numberofvertices, 1);
    md.hydrology.cavity_spacing = 2;
    md.hydrology.bump_height = 0.1*ones(md.mesh.numberofvertices, 1);
    md.hydrology.englacial_void_ratio = 1e-4;
    md.hydrology.omega = 1/2000;

    md.hydrology.sheet_alpha = 5./4.;
    md.hydrology.sheet_beta = 3./2.;

    % Allow channels and set channel conductivity
    md.hydrology.ischannels = 1;
    md.hydrology.channel_conductivity = 0.1;
    md.hydrology.channel_alpha = 1.25;
    md.hydrology.channel_beta = 1.5;

    % INITIAL CONDITIONS

    % Water layer thickness = 10% of bed bump height
    md.initialization.watercolumn = 0.1*md.hydrology.bump_height.*ones(md.mesh.numberofvertices, 1);
    md.initialization.channelarea = 0*ones(md.mesh.numberofedges, 1);

    % Set initial pressure equal to overburden
    phi_bed = md.constants.g*md.materials.rho_freshwater*md.geometry.base;
    p_ice = md.constants.g*md.materials.rho_ice*md.geometry.thickness;
    md.initialization.hydraulic_potential = phi_bed + p_ice;


    % BOUNDARY CONDITIONS
    % Set pressure=0 at terminus
    md.hydrology.spcphi = NaN(md.mesh.numberofvertices,1);
    pos=find(md.mesh.vertexonboundary & md.mesh.x==min(md.mesh.x));
    md.hydrology.spcphi(pos) = phi_bed(pos);

    % Specify no-flux Type 2 boundary conditions on all edges (except
    % the Type 1 condition set at the outflow above)
    md.hydrology.neumannflux=zeros(md.mesh.numberofelements,1);

    % FORCING
    md.hydrology.melt_flag = 1;
    md.basalforcings.groundedice_melting_rate = 0.05*ones(md.mesh.numberofvertices, 1);
    md.basalforcings.geothermalflux = 0;

    % Moulin inputs
    md.hydrology.moulin_input = moulin_inputs(dmesh.tri, 50, 0);
    save('MassCon.para.mat', 'md')
end 

if any(steps==3) 
    disp('	Step 3: Solve!');
    md=loadmodel('MassCon.para.mat');

    % Solve just hydrology
    md.transient=deactivateall(md.transient);
    md.transient.ishydrology=1;

    % Specify that you want to run the model on your current computer
    md.cluster=generic('np',4);

    % Define the time stepping scheme
    md.timestepping=timesteppingadaptive();
    md.timestepping.time_step_min=86400/md.constants.yts;
    md.timestepping.cfl_coefficient = 0.9;  % Must be <1 for stability
    md.timestepping.final_time = 10;
    md.settings.output_frequency = 30;

    md.initialization.vel = zeros(md.mesh.numberofvertices, 1) + 30;
    md.initialization.vx = zeros(md.mesh.numberofvertices, 1) - 30;
    md.initialization.vy = zeros(md.mesh.numberofvertices, 1) + 0;
    md.miscellaneous.name = 'MassCon';
    md = setmask(md,'','');

    % NUMERICAL PARAMETERS
    md.stressbalance.restol = 1e-3;
    md.stressbalance.reltol = nan;
    md.stressbalance.abstol = nan;
    md.stressbalance.maxiter = 100;

    md.verbose.solution=1;
    md=solve(md,'Transient');

    save('MassCon_transition_new.mat', 'md')
end 
