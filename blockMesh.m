% // Mesh generated on 08-Mar-2021
% // Local directory: /Users/fbisetti/teaching/ASE347/SP-2021/OF-assignments/of2/simulations/mesh
% // Username: fbisetti
% 
% 
inlet_height = 1;
theta = pi/6;
free_stream_length= 0.6;
tunnel_height = 0.6;
total_length = 3;

nx = 24;
ny = 8;
nz = 1;

vertecies(inlet_height,free_stream_length,tunnel_height, total_length, theta)
block(nx,ny,nz)
% blocks = block(nxCylinder,nyCylinder, nz, nHorizontalWake,nVertical,nHorizontalFore)
% [arcs] = arc(blocks,R)
% 
% inlets = inlet(blocks,v)
% outlets = outlet(blocks,v)
% tops = top (blocks,v)
% bottoms = bottom(blocks,v)
% cylinders = cylinder(blocks,v)
% 
% publish(blocks,v,arcs,inlets,outlets,tops,bottoms,cylinders)

% Now to pring out the file
function  [v] = vertecies(inlet_height,free_stream_length,tunnel_height, total_length, theta)
% setting the coordinate system to be (0,0,0) at the bottom left corner
% 
v = zeros(16,4);
% v(i,1) = x global coordinate
% v(i,2) = y global coordinate
% v(i,3) = z global coordinate
% v(i,4) = openFoam index 

% vertecies at inlet and freestream
for i = 1:2
    v(i,1) = (i-1)*free_stream_length;
    v(i+3,1) =  v(i,1);
    v(i+3,2) = inlet_height;
end

% vertecies at ramp
ramp_height = (inlet_height-tunnel_height)/2;
for i = 3
    v(i,1) =  free_stream_length+ramp_height/tan(theta);
    v(i,2) = ramp_height;
    v(i+3,1) =  v(i,1);
    v(i+3,2) = inlet_height-ramp_height;
end

% vertecies at outlet 
v(7,1) = total_length;
v(7,2) = ramp_height;

v(8,1) = total_length;
v(8,2) = inlet_height-ramp_height;


% Set vertecies for z axis 
v(9:16,1:2) = v(1:8,1:2);
for  i = 1:32
    v(i,3) = -0.05;
    v(i+32,3) = 0.05;    
end

% openFoam Indexing
for i = 1:16
    v(i,4) = i-1;
end

vx = v(:,1);
vy = v(:,2);
plot(vx,vy,'o');

end



function blocks = block(nx,ny,nz)
% nx number of cells in global x direction
% ny number of cells in global y direction
% nz number of cells in global z direction

blocks = zeros(3,15);
% blocks(i,1:8) vertex list
% blocks(i,9:11) nr of ceills in each local direction
% blocks(i,12:14) expansion ratio in each local direction
%  blocks(i,15) openFoam index

% inlet and ramp block
for i = 1:2
    % vertex list
    blocks(i,1) = i-1;
    blocks(i,2) = i;
    blocks(i,3) = blocks(i,2)+3;
    blocks(i,4) = blocks(i,1)+3;
    blocks(i,5) = blocks(i,1)+8;
    blocks(i,6) = blocks(i,2)+8;
    blocks(i,7) = blocks(i,3)+8;
    blocks(i,8) = blocks(i,4)+8;


end

% tunnel block 
% vertex list
blocks(3,1) = 2;
blocks(3,2) = 6;
blocks(3,3) = 7;
blocks(3,4) = 5;
blocks(3,5) = blocks(3,1)+8;
blocks(3,6) = blocks(3,2)+8;
blocks(3,7) = blocks(3,3)+8;
blocks(3,8) = blocks(3,4)+8;


% local cell number and simpleGrading 
blocks(1,9:11) = [nx,ny,nz]; 
blocks(2,9:11) = [nx,ny,nz]; 
blocks(3,9:11) = [4*nx,ny,nz];

blocks(1:3,12:14) = ones(3);

% OpenFoam block Index
for i = 1:3
    blocks(i,15) = i-1;
end
end


% 
% function [inlets] = inlet(blocks,v)
% inlets = [];
% left_end = min(v(:,1));
% for i = 1:(length(blocks))
%     if (v(blocks(i,1)+1,1) == left_end) & (v(blocks(i,4)+1,1) == left_end)
%         inlets = [inlets; blocks(i,1) blocks(i,4)+32 blocks(i,1)+32 blocks(i,4)];
%     end
% end
% end
% % function [outlets] = outlet(blocks,v)
% outlets = [];
% right_end = max(v(:,1));
% for i = 1:(length(blocks))
%     if (v(blocks(i,2)+1,1) == right_end) & (v(blocks(i,3)+1,1) == right_end)
%         outlets = [outlets; blocks(i,2) blocks(i,3)+32 blocks(i,2)+32 blocks(i,3)];
%     end
% end
% end

% function [tops] = top(blocks,v)
% 
% tops = [];
% top_end = max(v(:,2));
% for i = 1:(length(blocks))
%     if (v(blocks(i,3)+1,2) == top_end) & (v(blocks(i,4)+1,2) == top_end)
%         tops = [tops; blocks(i,3) blocks(i,4)+32 blocks(i,3)+32 blocks(i,4)];
%     end
% end
% end
% function [bottoms] = bottom(blocks,v)
% bottoms = [];
% bottom_end = min(v(:,2));
% for i = 1:(length(blocks))
%     if (v(blocks(i,1)+1,2) == bottom_end) & (v(blocks(i,2)+1,2) == bottom_end)
%         bottoms = [bottoms; blocks(i,1) blocks(i,1)+32 blocks(i,2)+32 blocks(i,2)];
%     end
% end
% end
% function [cylinders] = cylinder(blocks,v)
% cylinders = [];
% rad = 0.5;
% for i = 1:(length(blocks))
%     if (norm(v(blocks(i,1)+1,1:2),2) - rad < 0.0001) & (norm(v(blocks(i,4)+1,1:2),2) - rad < 0.0001)
%         cylinders = [cylinders; blocks(i,1) blocks(i,1)+8 blocks(i,2)+8 blocks(i,4)];
%     end
% end
% end
% 
% function publish(blocks,v,arcs,inlets,outlets,tops,bottoms,cylinders)
% 
% fid = fopen('test.txt','wt');
% % write header
% fprintf(fid,'// Mesh generated on 08-Mar-2021\n');
% fprintf(fid,'// Local directory: /Users/fbisetti/teaching/ASE347/SP-2021/OF-assignments/of2/simulations/mesh\n');
% fprintf(fid,'// Username: fbisetti\n\n\n');
% 
% fprintf(fid,'// Lf (fore)  =  4.00000e+00\n');
% fprintf(fid,'// Lw (wake)  =  6.00000e+00\n');
% fprintf(fid,'// R  (outer) =  1.00000e+00\n');
% fprintf(fid,'// H  (top/bottom) = 4.00000\n\n');
% 
% fprintf(fid,'FoamFile\n');
% fprintf(fid,'{\n\tversion  2.0;\n\tformat   ascii;\n\tclass    dictionary;\n\tobject   blockMeshDict;\n}\n\n');
% 
% fprintf(fid,'convertToMeters 1.0;\n\n');
% 
% % write verticies
% fprintf(fid,'verticies\n(\n');
% fprintf(fid,'\t( %+e %+e %+e ) // %1.0f\n',v(:,1:4)');
% fprintf(fid,');\n\n');
% 
% % write blocks
% fprintf(fid,'blocks\n(\n\n');
% fprintf(fid,'\t// block %1.0f \n\thex (%2.0f %2.0f %2.0f %2.0f %2.0f %2.0f %2.0f %2.0f) ( %2.0f %2.0f %2.0f) simpleGrading ( %+e %+e %3.1f)\n',[blocks(:,end) blocks(:,1:(end-1))]');
% fprintf(fid,'\n);\n\n');
% 
% % write arcs
% fprintf(fid,'edges\n(\n\n');
% fprintf(fid,'arc %2.0f %2.0f ( %+e %+e %+e )\n',arcs');
% fprintf(fid,'\n\n);\n\n');
% 
% % write boundaries
% fprintf(fid,'boundary\n(\n\n');
% 
% fprintf(fid,'\tinlet\n\t{\n\t\ttype patch;\n\t\tfaces\n\t\t(\n');
% fprintf(fid,'\t\t\t(%2.0f %2.0f %2.0f %2.0f)\n',inlets);
% fprintf(fid,'\t\t);\n\t}\n\n');
% 
% fprintf(fid,'\toutlet\n\t{\n\t\ttype patch;\n\t\tfaces\n\t\t(\n');
% fprintf(fid,'\t\t\t(%2.0f %2.0f %2.0f %2.0f)\n',outlets);
% fprintf(fid,'\t\t);\n\t}\n\n');
% 
% fprintf(fid,'\tcylinder\n\t{\n\t\ttype wall;\n\t\tfaces\n\t\t(\n');
% fprintf(fid,'\t\t\t(%2.0f %2.0f %2.0f %2.0f)\n',cylinders);
% fprintf(fid,'\t\t);\n\t}\n\n');
% 
% fprintf(fid,'\ttop\n\t{\n\t\ttype symmetryPlane;\n\t\tfaces\n\t\t(\n');
% fprintf(fid,'\t\t\t(%2.0f %2.0f %2.0f %2.0f)\n',tops);
% fprintf(fid,'\t\t);\n\t}\n\n');
% 
% fprintf(fid,'\tbottom\n\t{\n\t\ttype symmetryPlane;\n\t\tfaces\n\t\t(\n');
% fprintf(fid,'\t\t\t(%2.0f %2.0f %2.0f %2.0f)\n',bottoms);
% fprintf(fid,'\t\t);\n\t}\n\n');
% 
% fprintf(fid,');');
% 
% fclose(fid);
% 
% end