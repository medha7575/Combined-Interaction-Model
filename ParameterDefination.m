%% This files stores all the parameter values that are used in the respective simulation. 
%% It also creates the folder in which the data of the simulation will be stored at run time

cellCenters= [500,500; 300,500; 700,500; 400,330; 600,330; 400,670; 600,670; 200,330; 200, 670; 500,840; 800, 670; 800, 330; 500,160];



%% Following are the variables for easily assigning to array while initializing
TOTAL_NODES_AT_LOCATION = 3;
NODE_ID_AT_LOCATION=4;

ROW=1;
COL=2;
MYO_CONC=3;
NO_OF_NEIGHBOR=4;
MOVEORNOT=5;
COMBOVALUE=6;
MEAN_LENOFCONNEC=7;

BRANCHORNOT=4;
ASSCNPF=3;
THETA=5;
LEN=6;
MCSATBIRTH=7;
AGE=8;
%% Following are the parameters used 

% Lattice defination
bound=1000;
rows=bound;
cols=bound;

Combo = zeros(bound,bound); % Used for visualizing the lattice
npfChart= zeros(bound,bound); % Specifies where the nucleation promoters are present
Myo=zeros(bound,bound); % Used for visualizing the myosin activity field

Mc_NodeNo = ones(rows,cols).*-1; % Main array holding myosin material features
Ac_NodeNo = ones(rows,cols).*-1; % Main array holding actin material features


CapR1=20; % Inner radius of actin initialization
CapR=50;% Outer radius of actin initialization
MyoR=400;  % Whole myosin distribution area radius


Vnpf=50; %Value assigned to nucleation point in npfChart, Used for initialization
Ac= 150 ; % Actin value used for visualization
Mc= 350; % Active Myosin value for visualization
Dm=200; % Dead Myosin value for visualization


Den_c1=0.1;% initial npf density used for inner radius
Den_c=0.01; % initial npf density used for outer radius
Den_m=0.2; % myosin density

%% Other actin material parameters

AgeBr=9;  % age of branching
Knpf=0.5; % 'Kpoly' or Actin polymerization coefficient 
Lth=20; % Max length
AgeTh=20; % Age of Depolymerization 
P_del=0.7; %probability of deletion

%% Other actomyosin material parameters
Myo_conc=0.1; %Myosin activity per node/concnetration ('M'in paper)
Dis_Thres_Neigh1=5;%500nm distance Myosin Connection Search Distance ('Dthres' in paper)
MaxNeighAllowed=6; % Maximum number of neighbor allowed
Max_myosin_nodes_per_pixels =20;
K=0.5; %Coefficient of myosin pulling
KSp=0.05; %  Actin filamnet spring constant connecting myosin nodes


%% Intermaterial parameters
SR=5; %myosin and actin search radius for each other or Inter-material SR
Aura= 2; % Inter-material Volume Exclusion Threshold ('Vex' )



%% Defining folder and printing initial image
ver=0;  % used to create a new name for the file. Helpful in iterations
folder=strcat('FxHYwithCapDen_c1_',num2str(Den_c1),'Den_c_',num2str(Den_c),'Den_m_',num2str(Den_m),'Myo_conc',num2str(Myo_conc)); %this is where file will save
foldername=folder;
if exist(foldername,'dir')~=7
    mkdir(foldername);
else
    while  exist(foldername, 'dir')==7
        ver= ver+1;
        foldername= strcat(folder,num2str(ver));
    end
    mkdir(foldername)
end
absoluteFolderPath = foldername;

% % Setting directions used to initialize actin from npf. Used to direct
% arp 2/3 actin in different directions
k=1;
for i=-1:1
    for j=-1:1
        if i==0 & j==0
            continue
        else
            directions(k,1)= i; %r
            directions(k,2)= j; %c
            directions(k,3)= rad2deg(cart2pol(j,-i));
            k=k+1;
        end
    end
end
