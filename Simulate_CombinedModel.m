%% SIMULATION STARTS HERE
% This file is used to run the simulation after the initialization is
% completed. It calls module 'Force_Calculation' to first calculate forces on inidividual
% nodes of both material and then moves actin nodes by calling module
% 'MoveActin' and myosin nodes by calling module "MoveMysoin". 
% These modules are called randomly, so that both the materials have a random chance to
% get the first move.
% Thereafter, the code deletes filaments that have aged. It also
% removes wall to test affect of wall removal. Finally it generates output
% in form of images to visualize and study material behavior in the given
% setup.

NodeDensity = zeros(4000,2);
for mcs=1:2500
     mcs
   
   %%% initializing the arrays that are used in simulation
    NewActin=0;
    
    
    F_AcNet=zeros(Ac_NodeCount,2);
    
    Rho_myo=zeros(Ac_NodeCount,2);
    F_MyoNet=zeros(Mc_NodeCount,2); % total(spring + myosin) force on each node
    F_am=zeros(Mc_NodeCount,2);  % net Actomyosin force
    F_r_Nb= zeros(Mc_NodeCount,MaxNeighAllowed); %  Actomyosin force in x direction by ith neighbor
    F_c_Nb= zeros(Mc_NodeCount,MaxNeighAllowed);  %  actomyosin force in y direction by ith neighbor
    F_r_NbSp= zeros(Mc_NodeCount,MaxNeighAllowed); % Spring force in x direction by ith neighbor
    F_c_NbSp= zeros(Mc_NodeCount,MaxNeighAllowed);  % Spring force in y direction by ith neighbor
    F_Spring=zeros(Mc_NodeCount,2);     % Net PRing force by all nbeighbors
    F_myoT=zeros(Mc_NodeCount,2);
   
        %% % Forces calculation --- call associated module
  
       Force_MyosinThenActinNEW;
  
      
   Combo = zeros(bound,bound);
    R=rand(1);
	%% %Randomly chose actin or mysoin material for moving nodes
    
   if R>0.5
       MoveMyosin;
       MoveActin;
   else
       MoveActin;
       MoveMyosin;
   end
    
   
    %%%%%%%%%%%%%% For Deleting aged Filamnets  %%%%%%%%%%%%%%%%%%%%%%%
    kount=0;
    for NodeID= 1:Ac_NodeCount
        if NodeID>Ac_NodeCount
            break
        end
        if Ac_Node(NodeID,AGE) > (AgeTh)  %% instead of length compare age of actin
            Rno= rand(1);
            if Rno< P_del
                %             disp('count deleting nodes')
                kount= kount+1;
                NodeToDel(kount)= NodeID;
                
            end
        end
    end
        
for del = kount:-1:1  %% deleting in reverse order to avoid shifting overhead of nodeID
    NodeID= NodeToDel(del);
    Combo(round(Ac_Node(NodeID,ROW)),round(Ac_Node(NodeID,COL)))=0;
    %     disp('Deleteing Now')
    Ac_Node(NodeID,:)=[];
    Ac_NodeCount=Ac_NodeCount-1;
    
end
% % Reassigning Ac_NodeNo
Ac_NodeNo=ones(rows,cols).*-1;
for NodeID= 1:Ac_NodeCount
    Ac_NodeNo(round(Ac_Node(NodeID,ROW)),round(Ac_Node(NodeID,COL)))=NodeID;
end
    
    %%%%%%%%%%%%%%%%%%-----THE END--------%%%%%%%%%%%%%%%%%%%%%
    
     %%  VISUALIZATION and GENERATING IMAGES%%%%%
          
     if mcs==2
        close all
        cd (absoluteFolderPath)
        Actinfoldername=strcat('Combo');
        mkdir(Actinfoldername);
        cd (Actinfoldername)
        filename=strcat('Actin_',num2str(mcs));
        fig = figure;
        imagesc(Combo)
        axis square
        colormap('hot')
        
% % %              imshow(C, 'Colormap', jet(255))
        print(fig,filename,'-dpng');
        
        cd ..
        %% visualization of myosin concentration
        
        MyoCON=zeros(bound,bound);
        for i = 1:Mc_NodeCount
             roNo= round(Mc_Node(i ,ROW));
             coNo= round(Mc_Node(i ,COL));
             valu= Mc_Node(i ,MYO_CONC);
             MyoCON(roNo,coNo)=valu;
        end
        Myofoldername=strcat('MyoCon');
        mkdir(Myofoldername);
        cd (Myofoldername)
        
        filename=strcat('MyoCON',num2str(mcs));
        fig = figure;
        imagesc(MyoCON)
	    colormap('gray')
        print(fig,filename,'-dpng');
                   
        
        
        cd ..
        cd ..
     end
     
%     if mod(mcs,100)==0
%         
%         tbl=tabulate( Mc_Node(:,NO_OF_NEIGHBOR));
%         cd (absoluteFolderPath)
%         filename=strcat('FreqN_',num2str(mcs));
%         fig = figure;
%         bar(tbl(:,1),tbl(:,2))
%          print(fig,filename,'-dpng');
%         cd ..
%     end
    
    if mod(mcs,50)==0
        close all
        cd (absoluteFolderPath)
        Actinfoldername=strcat('Combo');
        mkdir(Actinfoldername);
        cd (Actinfoldername)
        filename=strcat('Actin_',num2str(mcs));
        fig = figure;
        imagesc(Combo)
        axis square
        colormap('hot')
        
        %      imshow(C, 'Colormap', jet(255))
        print(fig,filename,'-dpng');
        
        cd ..
                

        MyoCON=zeros(bound,bound);
        for i = 1:Mc_NodeCount
             roNo= round(Mc_Node(i ,ROW));
             coNo= round(Mc_Node(i ,COL));
             valu= Mc_Node(i ,MYO_CONC);
             MyoCON(roNo,coNo)=valu;
        end
        cd (Myofoldername)
        filename=strcat('MyoCON',num2str(mcs));
        fig = figure;
        imagesc(MyoCON)
        colormap('gray')
        print(fig,filename,'-dpng');
                   
        
        cd ..
        cd ..
        
%         
    end
    %
    M=Ac_NodeNo>0;
    nnz(M);
 ActinLenCounter(mcs)=length(find(Ac_Node(:,6)<=AgeBr));
 ActinCounter(mcs)=Ac_NodeCount;
 NewBranches(mcs)=NewActin;
 DelActin(mcs)=kount;
 
 
 
end
   cd (absoluteFolderPath)
    filename=strcat('ActotalPlot_',num2str(mcs));
    fig = figure;
    plot(ActinCounter)
     print(fig,filename,'-dpng');
     
    filename=strcat('LenPlot_',num2str(mcs));
    fig = figure;
    plot(ActinLenCounter)
     print(fig,filename,'-dpng');
     
         filename=strcat('NewActinPlot_',num2str(mcs));
    fig = figure;
    plot(NewBranches)
     print(fig,filename,'-dpng');
     
         filename=strcat('DelActinPlot_',num2str(mcs));
    fig = figure;
    plot(DelActin)
     print(fig,filename,'-dpng');
    cd ..

