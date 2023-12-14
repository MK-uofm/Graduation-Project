gm =importGeometry("./Large_free.STL");

% model = createpde(1);
% b=importGeometry(model,'./Small_free.STL');
% g=model.Geometry;

model.Geometry = gm;
e = nearestEdge(gm, [0 , 0, 0 ]);
facesAttachedToEdges(gm, e);
% vertices = gm.vertexCoordinates(1: gm.NumVertices);
% pdegplot(gm,FaceLabels="on", FaceAlpha=0.3, EdgeLabels="on");
[Faces, Vertices] = gm.allDisplayFaces();

n = length(Faces);
Edgelist = zeros(n*3,2);
curr = 1;
for i=1:n
    Edgelist(curr,1) = min(Faces(i,1),Faces(i,2));
    Edgelist(curr,2) = max(Faces(i,1),Faces(i,2));
    Edgelist(curr+1,1) = min(Faces(i,2),Faces(i,3));
    Edgelist(curr+1,2) = max(Faces(i,2),Faces(i,3));
    Edgelist(curr+2,1) = min(Faces(i,1),Faces(i,3));
    Edgelist(curr+2,2) = max(Faces(i,1),Faces(i,3));
    curr = curr + 3;
end
numV = length(Vertices);
for i=1:numV
    d(i) = sqrt(Vertices(i,1)^2+Vertices(i,2)^2);
end
EList = unique(Edgelist,'rows');
EdgeProperties = zeros(length(EList),5);
for i = 1:length(EList) 
    EdgeProperties(i,1) = d(EList(i,1));
    EdgeProperties(i,2) = d(EList(i,2));
    EdgeProperties(i,3) = Vertices(EList(i,1),3);
    EdgeProperties(i,4) = Vertices(EList(i,2),3);
    EdgeProperties(i,5) = sqrt((Vertices(EList(i,1),1)-Vertices(EList(i,2),1))^2+(Vertices(EList(i,1),2)-Vertices(EList(i,2),2))^2);
end

isSameDistance = (abs(EdgeProperties(:,1)-EdgeProperties(:,2))<0.00001);
isSameHeight = (abs(EdgeProperties(:,3)-EdgeProperties(:,4))<0.00001);
isSameCircle = and(isSameHeight,isSameDistance);
indSameCircle = find(isSameCircle);
isInner = and(isSameDistance,(EdgeProperties(:,1)<21));
isOuter = and(isSameDistance,~isInner);
isBottom = and(isSameHeight,(EdgeProperties(:,3)<1));
isTop = and(isSameHeight,~isBottom);
isInnerBottom = and(isBottom, isInner);
indInnerBottom = find(isInnerBottom);
n = length(indInnerBottom);
theta = pi()/n;
R = mean(EdgeProperties(indInnerBottom(:),1));
Atri = R^2*cos(theta)*sin(theta);
Ahole = n * Atri;
Acircle = pi()*R^2;
h = R*(1-cos(theta));

disp("Num Edges: ");
disp(length(isTop));
disp("Num Vertices: ");
disp(numV);
disp("Num Faces: ");
disp(length(Faces));

Dpoly = 2*(R-h);
d_hole_min = 40;
d_hole_max = 40.025;
d_shaft_min = 39.984;
d_shaft_max = 40.000;

adjusted_r_hole_min = (d_hole_min+2*h)/2;
adjusted_r_hole_max = (d_hole_max+2*h)/2;
adjusted_r_shaft_min = (d_shaft_min+2*h)/2;
adjusted_r_shaft_max = (d_shaft_max+2*h)/2;


fprintf('Type of fit: Locational Clearance\n')
fprintf('Diameter of "polygon":');
disp(Dpoly)
fprintf('Hole diameter min:');
disp(2*adjusted_r_hole_min);
fprintf('Hole diameter max:');
disp(2*adjusted_r_hole_max);

fprintf('Shaft diameter min:');
disp(2*adjusted_r_shaft_min);
fprintf('Shaft diameter max:');
disp(2*adjusted_r_shaft_max);



