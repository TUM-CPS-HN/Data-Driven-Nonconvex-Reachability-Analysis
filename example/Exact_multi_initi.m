function Z = Copy_2_of_exact_multi(Z1,Z2,system_info,varargin)
% exactAnd - overloads & operator, computes the multiplication of constraint matrix zonotope and constraint poly zonotopes
%
% Syntax:
%    Z = exactAnd(Z1,Z2)
%
% Inputs:
%    Z1 - constraint zonotope
%    Z2 - constraint poly zonotope,
%
% Outputs:
%    Z - new constraint poly zonotope
%
% Example:
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Zhen Zhang, Amr Alanwar
% Written:       11-Mar-2025
% Last update:
%
%
% Last revision: ---

%------------- BEGIN CODE --------------
nx=system_info.nx;
nu=system_info.nu;

Z1C = reshape(Z1.c, nx, nx+nu);
Z1G = {};
for i=1:size(Z1.G,2)
    Z1G{i} = reshape(Z1.G(:,i), nx, nx+nu);
end

    Z1A = reshape(Z1.A, [], 1);

Z1B = reshape(Z1.b, [], 1);

if(~isempty(Z1C) && ~isempty(Z2.c))
    newcen = Z1C * Z2.c;
end




zeroVec = zeros(size(Z1C,1));
% indexI =1;
newE =[];
newA=[];
newb=[];
newEC=[];
index =1;
newGen ={};

% bring the exponent matrices to a common representation
[idCom,E1Com,E2Com] = mergeExpMatrix(Z1.id,Z2.id,Z1.E,Z2.E);
[idComC,EC1Com,EC2Com] = mergeExpMatrix(Z1.id,Z2.id,Z1.EC,Z2.EC);

[rE1,cE1] = size(E1Com);
[rE2,cE2] = size(E2Com);

if rE1 >rE2
    E2ComZeros = [E2Com;zeros(rE1-rE2,cE2)];
    E1ComZeros = E1Com;
else
    E1ComZeros = [E1Com;zeros(rE2-rE1,cE1)];
    E2ComZeros = E2Com;
end

E1Com = E1ComZeros;
E2Com = E2ComZeros;


[rEC1,cEC1] = size(EC1Com);
[rEC2,cEC2] = size(EC2Com);

if rEC1 >rEC2
    EC2ComZeros = [EC2Com;zeros(rEC1-rEC2,cEC2)];
    EC1ComZeros = EC1Com;
else
    EC1ComZeros = [EC1Com;zeros(rEC2-rEC1,cEC1)];
    EC2ComZeros = EC2Com;
end

EC1Com = EC1ComZeros;
EC2Com = EC2ComZeros;


%G1 * c2
if(~isempty(Z2.c) && ~isempty(Z1G))
    for i=1:size(Z1G,2)
        vecMultip =  Z1G{i} * Z2.c;
        if isequal(vecMultip,zeroVec)
            continue;
        end
        newGen{index} = vecMultip;
        newE = [newE E1Com(:,i)];
        index=index+1;
    end
end

%c1 * G2
if(~isempty(Z1C) && ~isempty(Z2.G))
    for i=1:size(Z2.G,2)
        vecMultip = Z1C * Z2.G(:,i);
        if isequal(vecMultip,zeroVec)
            continue;
        end
        newGen{index} = vecMultip;
        newE = [newE E2Com(:,i)];        
        index=index+1;
    end
end

%G1 * G2
if(~isempty(Z1G) && ~isempty(Z2.G))
    for i=1:length(Z1G)
        for k=1:size(Z2.G,2)
            vecMultip = (Z1G{i} * Z2.G(:,k));
            if isequal(vecMultip,zeroVec)
                continue;
            end
            newGen{index} = vecMultip;
            tempECom = E1Com(:,i) + E2Com(:,k);
            newE = [newE tempECom]; 
            index=index+1;
        end
    end
end

% A
tempnewA=[];

    newA = Z2.A                               ;



    newb = [reshape(Z1B, [], 1);
            Z2.b];

% EC
if(~isempty(Z1.EC) && ~isempty(1))
    newEC = [EC1Com EC2Com];
end



MatnewGen=[];
for i=1:size(newGen,2)
    MatnewGen=[MatnewGen newGen{i}];
end

Z = conPolyZono(newcen,MatnewGen,newE,newA,newb,newEC,[],idCom);


end

%------------- END OF CODE --------------