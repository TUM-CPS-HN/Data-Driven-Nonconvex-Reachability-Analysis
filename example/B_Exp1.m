% t_linearDT - computes the data driven reachable set of discrete time systems
% x(k+1) = Ax(k) + Bu(k) + w(k)
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: Amr Alanwar, Anne Koch, Frank Allgöwer, Karl Johansson "Data Driven Reachability Analysis Using Matrix Zonotopes"
%
%
%
% Author:       Zhen Zhang, Amr Alanwar
% Written:      13-March-2025
% Last update:  
% Last revision:---

%------------- BEGIN CODE --------------

rand('seed',1);

clear all
% close all
%% system dynamics
dim_x = 5;
A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B_ss = ones(5,1);
nx = size(A, 1);
nu = size(B_ss, 2);

C = [1,0,0,0,0];
D = 0;
% define continuous time system
sys_c = ss(A,B_ss,C,D);
% convert to discrete system
samplingtime = 0.05;
sys_d = c2d(sys_c,samplingtime);


%Number of trajectories
initpoints =2;

%Number of trajectories
steps1 = 25; % new data start point steps = 200 steps1 = 100
steps = 100;
totalsamples = initpoints*steps;
totalsamples1 = initpoints*steps1;
system_info = struct('nx',nx, 'nu',nu);

%% initial set and input
X0 = zonotope(ones(dim_x,1),0.1*diag(ones(dim_x,1)));
U = zonotope(10,0.25);
scale=1;
%noise zontope W
W = zonotope(zeros(dim_x,1),scale*0.005*ones(dim_x,1));

W_new = zonotope(zeros(dim_x,1),scale*0.005*ones(dim_x,1));

%Construct matrix zonotpe \mathcal{M}_w
index=1;
for i=1:size(W.generators,2)
    vec=W.Z(:,i+1);
    GW{index}= [ vec,zeros(dim_x,totalsamples1-1)];
    for j=1:totalsamples1-1
        GW{j+index}= [GW{index+j-1}(:,2:end) GW{index+j-1}(:,1)];
    end
    index = j+index+1;
end
Wmatzono= matZonotope(zeros(dim_x,totalsamples1),GW);


index=1;
for i=1:size(W_new.generators,2)
    vec=W_new.Z(:,i+1);
    GW_new{index}= [ vec,zeros(dim_x,(totalsamples - totalsamples1 -1))];
    for j=1:(totalsamples - totalsamples1 -1)
        GW_new{j+index}= [GW_new{index+j-1}(:,2:end) GW_new{index+j-1}(:,1)];
    end
    index = j+index+1;
end
Wmatzono_new= matZonotope(zeros(dim_x,( totalsamples - totalsamples1)),GW_new);

index=1;
for i=1:size(W_new.generators,2)
    vec=W_new.Z(:,i+1);
    GW_all{index}= [ vec,zeros(dim_x,(totalsamples -1))];
    for j=1:(totalsamples -1)
        GW_all{j+index}= [GW_all{index+j-1}(:,2:end) GW_all{index+j-1}(:,1)];
    end
    index = j+index+1;
end
Wmatzono_all= matZonotope(zeros(dim_x,( totalsamples)),GW_all);



% randomly choose constant inputs for each step / sampling time
for i=1:totalsamples
    u(i) = randPoint(U);
end


%simulate the system to get the data
x0 = X0.center;
x(:,1) = x0;
index=1;
for j=1:dim_x:initpoints*dim_x
    x(j:j+dim_x-1,1) = randPoint(X0);
    for i=1:steps
        utraj(j,i) = u(index);
        x(j:j+dim_x-1,i+1) = sys_d.A*x(j:j+dim_x-1,i) + sys_d.B*u(index) + randPoint(W);      
        index=index+1;
    end
end

% new %

% concatenate the data trajectories 
index_0 =1;
index_1 =1;
for j=1:dim_x:initpoints*dim_x
    for i=2:steps+1
        x_meas_vec_1(:,index_1) = x(j:j+dim_x-1,i);
        index_1 = index_1 +1;
    end
    for i=1:steps
        u_mean_vec_0(:,index_0) = utraj(j,i);
        x_meas_vec_0(:,index_0) = x(j:j+dim_x-1,i);
        index_0 = index_0 +1;
    end
end

U_full = u_mean_vec_0(:,1:totalsamples1); %same as u 
X_0T = x_meas_vec_0(:,1:totalsamples1);
X_1T = x_meas_vec_1(:,1:totalsamples1);

U_full_new = u_mean_vec_0(:,totalsamples1+1:totalsamples); %same as u 
X_0T_new = x_meas_vec_0(:,totalsamples1+1:totalsamples);
X_1T_new = x_meas_vec_1(:,totalsamples1+1:totalsamples);

U_full_all = u_mean_vec_0(:,1:totalsamples); %same as u 
X_0T_all = x_meas_vec_0(:,1:totalsamples);
X_1T_all = x_meas_vec_1(:,1:totalsamples);

X1W_cen =  X_1T - Wmatzono.center;
X1W = matZonotope(X1W_cen,Wmatzono.generator);

% set of A and B
AB = X1W  *pinv([X_0T;U_full]);

% validate that A and B are within AB
intAB11 = intervalMatrix(AB);
intAB1 = intAB11.int;
intAB1.sup >= [sys_d.A,sys_d.B]
intAB1.inf <= [sys_d.A,sys_d.B]


X1W_cen_all =  X_1T_all - Wmatzono_all.center;
X1W_all = matZonotope(X1W_cen_all,Wmatzono_all.generator);

% set of A and B
AB_all = X1W_all  *pinv([X_0T_all;U_full_all]);

% validate that A and B are within AB
intAB11_all = intervalMatrix(AB_all);
intAB1_all = intAB11_all.int;
intAB1_all.sup >= [sys_d.A,sys_d.B]
intAB1_all.inf <= [sys_d.A,sys_d.B]

C1=(X_1T - Wmatzono.center)*pinv([X_0T;U_full]);

for i=1:totalsamples1
    
    G1(:,:,i)=Wmatzono.generator(:,:,i)*pinv([X_0T;U_full]);

end

M1=matZonotope(C1,G1);
M1=zonotope(M1);
M1 = conZonotope(M1);

C2=(X_1T_new - Wmatzono_new.center)*pinv([X_0T_new;U_full_new]);

for i=1:totalsamples - totalsamples1
    
    G2(:,:,i)=Wmatzono_new.generator(:,:,i)*pinv([X_0T_new;U_full_new]);

end

M2=matZonotope(C2,G2);
M2=zonotope(M2);
M2 = conZonotope(M2);

MM3=ConMat_Intersection(M1,M2);
M3=MM3;
C3=MM3.c;
G3=MM3.G;
A3=MM3.A;
B3=MM3.b;
E3=eye(size(G3,2));
EC3=eye(size(E3,1));

polyM3 = conPolyZono(C3,G3,E3,A3,B3,EC3);

polyX = conPolyZono(X0);
polyU = conPolyZono(U);
polyW = conPolyZono(W);

%% compute next step sets from model / data

% set number of steps in analysis
totalsteps = 5;
X_model = cell(totalsteps+1,1);
X_data = cell(totalsteps+1,1);
X_exact = cell(totalsteps+1,1);
X_all = cell(totalsteps+1,1);
X_exact2 = cell(totalsteps+1,1);
% init sets for loop
X_model{1} = X0; X_data{1} = X0; X_exact{1} = polyX; X_all{1} = X0;

redOp=800;
for i=1:totalsteps
    
    % 1) model-based computation
    
    X_model{i+1,1} = sys_d.A * X_model{i} + sys_d.B * U+W;
    X_model{i+1,1}=reduce(X_model{i+1,1},'girard',redOp);
    % 2) Data Driven approach
   
    X_data{i+1,1} = AB * (cartProd(X_data{i},U)) +W;
    X_data{i+1,1}=reduce(X_data{i+1,1},'girard',redOp);

    X_all{i+1,1} = AB_all * (cartProd(X_all{i},U)) +W;
    X_all{i+1,1}=reduce(X_all{i+1,1},'girard',redOp);


    % 3) exact multiplication

    P = cartProd(X_exact{i},polyU);
    X_exact{i+1,1} = exactPlusWConst(Exact_Multi(polyM3,P,system_info),polyW);
    X_exact{i+1,1}=reduce(X_exact{i+1,1},'girard',redOp);

    % 

end



%% visualization

projectedDims = {[1 2],[3 4],[4 5]};
axx{1} = [0.75,1.5,0.5,4]; axx{2} = [0.75,3,0.8,2.2];axx{3} = [0.75,2.3,0.75,2.8];
index=1;
numberofplots = totalsteps;%length(X_model)
for plotRun=1:length(projectedDims)
    
    figure('Renderer', 'painters', 'Position', [10 10 700 900])
    
    % set axis
    % plotRun=3; 
    index=index+1;
    % plot initial set
    handleX0 = plot(X0,projectedDims{plotRun},'k-','LineWidth',2);
    hold on;
    
    % iSet=5;
    % plot reachable sets starting from index 2, since index 1 = X0
    
    % plot reachable sets from model
    for iSet=2:numberofplots
    handleModel=  plot(X_model{iSet},projectedDims{plotRun},'b','Filled',true,'FaceColor',[.8 .8 .8],'EdgeColor','b');
    end
    
    % plot reachable sets from data
    for iSet=2:numberofplots
    handleData=   plot(X_data{iSet},projectedDims{plotRun},'b');
    end

   

    % plot reachable sets from exact
    for iSet=2:numberofplots
        handleExact=   plot(X_exact{iSet},projectedDims{plotRun},'r');
    end

    for iSet=2:numberofplots
        handleAll=   plot(X_all{iSet},projectedDims{plotRun},'k');
    end

    % break;

    % label plot
    xlabel(['x_{',num2str(projectedDims{plotRun}(1)),'}']);
    ylabel(['x_{',num2str(projectedDims{plotRun}(2)),'}']);
    %axis(axx{plotRun});
    % skip warning for extra legend entries
    warOrig = warning; warning('off','all');
    legend([handleX0,handleModel,handleExact,handleAll,handleData],...
        {'Initial Set','Model Based','$\hat{\mathcal{R}}$ Algorithm 1', '$\hat{\mathcal{R}}_a$ (Alanwar et al.)','$\hat{\mathcal{R}}_1$ Initial Reachability'}, ...
        'Interpreter', 'latex', 'Location', 'best');
    warning(warOrig);
    ax = gca;
    ax.FontSize = 22;
    %set(gcf, 'Position',  [50, 50, 800, 400])
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

end

%------------- END OF CODE --------------

