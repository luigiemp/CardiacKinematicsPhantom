function PlotToVTK(xt, conn, ElType, NodeScalars, NodeScalarsName, CellScalars, CellScalarsName, OutputFileName, TimeStep)
% Print connectivity, invariants, and nodal positions to VTK

NodeNumber = size(xt,1);
ElementNumber = size(conn,1);
NodePerElement = size(conn,2);

% Printing vtk file
fileID = fopen([OutputFileName,num2str(TimeStep),'.vtk'],'w');
fprintf(fileID, '%s\n', '# vtk DataFile Version 4.0');
fprintf(fileID, '%s\n', 'vtk output');
fprintf(fileID, '%s\n', 'ASCII');
fprintf(fileID, '%s\n', 'DATASET UNSTRUCTURED_GRID');
fprintf(fileID, '%s ' , 'POINTS');
fprintf(fileID, '%d ' , NodeNumber);
fprintf(fileID, '%s\n', 'float');

for i = 1:NodeNumber
    fprintf(fileID, '%1.6f %1.6f %1.6f\n', xt(i,:));
end

fprintf(fileID, '\n%s ' , 'CELLS');
fprintf(fileID, '%d ' , ElementNumber);
fprintf(fileID, '%d\n', ElementNumber*(NodePerElement+1));



switch ElType
    case 11 % triangles + linear
        for i = 1:ElementNumber
            fprintf(fileID, '%d %d %d %d\n',[3, conn(i,1:3)]);
        end
        CellType = 5;
    case 12 % triangles + quadratic
        for i = 1:ElementNumber
            fprintf(fileID, '%d %d %d %d %d %d %d\n',[6, conn(i,1:6)]);
        end
        CellType = 22;
    case 21 % quadrilateral + linear
        for i = 1:ElementNumber
            fprintf(fileID, '%d %d %d %d %d\n',[4, conn(i,1:4)]);
        end
        CellType = 9;
    case 22 % quadrilateral + quadratic
        for i = 1:ElementNumber
            fprintf(fileID, '%d %d %d %d %d %d %d %d %d\n',[8, conn(i,1:8)]);
        end
        CellType = 23;
    case 31 % tetrahedral + linear
        for i = 1:ElementNumber
            fprintf(fileID, '%d %d %d %d %d\n',[4, conn(i,1:4)]);
        end
        CellType = 10;
    case 32 % hexahedral + linear
        for i = 1:ElementNumber
            fprintf(fileID, '%d %d %d %d %d %d %d %d %d\n',[8, conn(i,1:8)]);
        end
        CellType = 12;
    otherwise
        error('Function PlotToVTK: Connectivity type not implemented.');
end

fprintf(fileID, '\n%s ' , 'CELL_TYPES');
fprintf(fileID, '%d\n' , ElementNumber);

for e = 1:ElementNumber
    fprintf(fileID, '%d\n',CellType);
end



% Print node scalars
if (size(NodeScalarsName,1) > 0)
    fprintf(fileID, '\n%s ' , 'POINT_DATA ');
    fprintf(fileID, '%d\n' , NodeNumber);

    for sc = 1:size(NodeScalars, 2)
        fprintf(fileID, '%s\n', ['SCALARS ',NodeScalarsName(sc,:),' double']);
        fprintf(fileID, '%s\n', 'LOOKUP_TABLE default');
        for n = 1:NodeNumber
            fprintf(fileID, '%1.6f\n', NodeScalars(n,sc) );
        end
    end
end



% Print cell scalars
if (size(CellScalarsName,1) > 0)
    fprintf(fileID, '\n%s ' , 'CELL_DATA ');
    fprintf(fileID, '%d\n' , ElementNumber);

    for sc = 1:size(CellScalars,2)
        fprintf(fileID, '%s\n', ['SCALARS ',CellScalarsName(sc,:),' double']);
        fprintf(fileID, '%s\n', 'LOOKUP_TABLE default');
        for e = 1:ElementNumber
            fprintf(fileID, '%1.6f\n', CellScalars(e,sc) );
        end
    end
end



fclose(fileID);



end % end of PlotToVTK function













