function [] = VectorsToVTK(Nodes, Vectors, VectorName, Scalars, ScalarsName, OutputFileName, time)
                      
% File name
fileID = fopen([OutputFileName,'_',num2str(time),'.vtk'],'w');

NumNodes = size(Nodes,1);

% Header file
fprintf(fileID, '%s\n', '# vtk DataFile Version 4.0');
fprintf(fileID, '%s\n', 'vtk output');
fprintf(fileID, '%s\n', 'ASCII');
fprintf(fileID, '%s\n', 'DATASET POLYDATA');
fprintf(fileID, '%s ' , 'POINTS');
fprintf(fileID, '%d ' , NumNodes*2);
fprintf(fileID, '%s\n', 'float');

% Nodal positions
if (size(Nodes,2) == 2)
    Nodes(:,3) = 0;
end

for i = 1:NumNodes
    fprintf(fileID, '%1.6f %1.6f %1.6f \n', Nodes(i,:));
    fprintf(fileID, '%1.6f %1.6f %1.6f \n', Nodes(i,:)); % Twice so Vector is centered at nodes
end

% Scalar values 
if (size(Scalars, 1) > 0)
    fprintf(fileID, '\n%s ' , 'POINT_DATA');
    fprintf(fileID, '%d \n' , NumNodes*2 );
    for sc = 1:size(Scalars,2)
        fprintf(fileID, '%s \n' , ['SCALARS   ',ScalarsName(sc,:),'   double']);
        fprintf(fileID, '%s\n', 'LOOKUP_TABLE default');
        for i = 1:size(Scalars,1)
            fprintf(fileID, '%1.6f \n', Scalars(i, sc));
            fprintf(fileID, '%1.6f \n', Scalars(i, sc));
        end
    end
end

if (size(Scalars, 1) == 0)
    fprintf(fileID, '\n%s ' , 'POINT_DATA');
    fprintf(fileID, '%d \n' , NumNodes*2);
end
fprintf(fileID, '%s \n' , ['VECTORS   ',VectorName,'   double']);

% Vector directions
for i = 1:NumNodes
    fprintf( fileID, '%1.6f %1.6f %1.6f \n',  0.5*Vectors(i,:) );
    fprintf( fileID, '%1.6f %1.6f %1.6f \n', -0.5*Vectors(i,:) );
end

fclose(fileID);

end


