function res = ConMat_Intersection(cZ,S)

%    % constrained zonotopes
%    Z = [0 3 0 1;0 0 2 1];
%    A = [1 0 1]; b = 1;
%    cZ1 = conZonotope(Z,A,b);
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1]; b = 1;
%    cZ2 = conZonotope(Z,A,b);
%
%    % halfspace and constrained hyperplane
%    P1 = polytope([1,-2],1);
%    P2 = polytope([-2 -0.5;1 0],[-4.25;2.5],[1,-2],1);
%
%    % compute intersection
%    res1 = cZ1 & cZ2;
%    res2 = cZ2 & P1;
%    res3 = cZ1 & P2;
%
%    % visualization
%    figure; hold on;
%    plot(cZ1,[1,2],'r');
%    plot(cZ2,[1,2],'b');
%    plot(res1,[1,2],'FaceColor','g');
%    title('Constrained zonotope');
%
%    figure; hold on; xlim([-4,4]); ylim([-4,4]);
%    plot(P1,[1,2],'r','FaceAlpha',0.5);
%    plot(res2,[1,2],'FaceColor','g');
%    plot(cZ2,[1,2],'b');
%    title('halfspace');
%
%    figure; hold on; xlim([0,4]); ylim([-3,4]);
%    plot(P2,[1,2],'g');
%    plot(cZ1,[1,2],'r');
%    plot(res3,[1,2],'b','LineWidth',2);
%    title('Constrained hyperplane'); 


% Calculate intersection according to equation (13) at Proposition 1 in
% reference paper [1]
Z = [cZ.c, cZ.G, zeros(size(S.G))];
A = blkdiag(cZ.A,S.A);
A = [A; cZ.G, -S.G];
b = [cZ.b; S.b; S.c - cZ.c];

res = conZonotope(Z,A,b);

% delete all zero constraints and generators
% res = compact_(res,'zeros',eps);

end

% ------------------------------ END OF CODE ------------------------------
