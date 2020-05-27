close all
clear

%read building summary spreadsheet
bldg= readmatrix('bldg_table_removed small area.xlsx');

%read original building spreadpsheet
original_bldg= readmatrix('bldg_removed no wall.xlsx');

for bldg_num= 1:1
    %extract information needed to create .key file
    bldg_id= bldg(bldg_num,21);
    
    %create file    
    key_file= fopen(sprintf('Bldg%d_Keyword_File.key',bldg_id),'w');
    fprintf(key_file,'*KEYWORD\n');

    %read processed wall data spreadsheet
    wall= readmatrix(sprintf('Bldg Plan Jan 4/Wall Table/wall_table%d.xlsx',bldg_id));
    wall(isnan(wall))=0;
    
    %read wall orientation from Excel
    wall_orientation= readcell(sprintf('Bldg Plan Jan 4/Wall Table/wall_table%d.xlsx',bldg_id));
    wall_orientation= wall_orientation(2:end,8);
    
    %*MAT_ELASTIC
    
    %initialize matrix
    conc_mat= zeros(1,7); 
    
    %input values to matrix
    conc_mat(1,1)= 1;                                                       %material ID
    conc_mat(1,2)= 2400;                                                    %mass density
    %Young's modulus
    if isnan(original_bldg(bldg_num,48))==1
        conc_mat(1,3)= 4.5*10^9*sqrt(4000*0.00689476);
    else
        conc_mat(1,3)= 4.5*10^9*sqrt(original_bldg(bldg_num,48)*0.00689476);
    end
    conc_mat(1,4)= 0.17;                                                    %Poisson's ratio
    conc_mat(1,5)= 0;                                                       %axial damping factor
    conc_mat(1,6)= 0;                                                       %bending damping factor
    conc_mat(1,7)= conc_mat(1,3)/(3*(1-2*conc_mat(1,4)));                   %bulk modulus
    
    %print heading
    fprintf(key_file,'*MAT_ELASTIC\n$Material ID,Mass Density,Youngs Modulus,Poissons Ratio,Axial Damping Factor,Bending Damping Factor,Bulk Modulus\n');
    
    %print matrix
    for i=1:size(conc_mat,2)
        if i<size(conc_mat,2)
            fprintf(key_file, '%.4g,', conc_mat(1,i));
        elseif i== size(conc_mat,2)
            fprintf(key_file, '%.4g', conc_mat(1,i));
            fprintf(key_file, '\n');
        end
    end
    
    %*NODE    
    height= bldg(bldg_num,5);
    
    num_story= bldg(bldg_num,4);
    
    %specify default constraint conditions (see manual)
    node_tc= 0;
    node_rc= 0;
    
    node_data= zeros(size(wall,1)*(num_story+1),6);
    
    node_data(:,5)= node_tc;
    node_data(:,6)= node_rc;
    
    row= 1;
    
    for story=1:num_story+1
        node_data(row:row+size(wall,1)-1,2:3)= wall(:,2:3);
        node_data(row:row+size(wall,1)-1,1)= story*1000+[1:size(wall,1)]';
        node_data(row:row+size(wall,1)-1,4)= (story-1)*height;
        
        row= row+ size(wall,1);
    end
    
    fprintf(key_file,'*NODE\n$Node Label,X Coordinate,Y Coordinate,Z Coordinate,TC,RC\n');
    
    for i =1:size(node_data,1)
        for ii=1:size(node_data,2)
            if ii== size(node_data,2)
                fprintf(key_file, '%d', node_data(i,ii));
            elseif ii< size(node_data,2)
                fprintf(key_file, '%d,', node_data(i,ii));
            end
        end
        fprintf(key_file, '\n');
    end
    
    %node 3 for elements
    node3_data= zeros(size(wall,1),6);
    
    %place node 3 at 0 elevation
    node3_data(:,1)= [1:size(wall,1)]';
    node3_data(:,2)= wall(:,2);
    node3_data(:,3)= wall(:,2)+10000;
    node3_data(:,4)= 0;
    node3_data(:,5)= 7;
    node3_data(:,6)= 7;
    
    fprintf(key_file, '$Element Node 3\n');
    
    for i =1:size(node3_data,1)
        for ii=1:size(node3_data,2)
            if ii== size(node3_data,2)
                fprintf(key_file, '%d', node3_data(i,ii));
            elseif ii< size(node3_data,2)
                fprintf(key_file, '%d,', node3_data(i,ii));
            end
        end
        fprintf(key_file, '\n');
    end
    
    %mass node calculation
    %1.area 2.mass 3.x*mass 4.y*mass
    %all standard metric
    mass_calc= zeros(size(wall,1),4);
    
    story_height= bldg(bldg_num,5);
    
    mass_calc(:,1)= wall(:,4).*wall(:,5);
    mass_calc(:,2)= conc_mat(1,2)*story_height*mass_calc(:,1);
    mass_calc(:,3)= wall(:,2).* mass_calc(:,2);
    mass_calc(:,4)= wall(:,3).* mass_calc(:,2);
    
    overall_x_com= bldg(bldg_num,12);
    overall_y_com= bldg(bldg_num,13);
    
    seismic_mass= bldg(bldg_num,3)*(bldg(bldg_num,9)*conc_mat(1,2)+ 1200/9.81);
    
    %back calculate slab COM given available data
    slab_x_com= (overall_x_com*(seismic_mass+sum(mass_calc(:,2)))- sum(mass_calc(:,3)))/seismic_mass;
    slab_y_com= (overall_y_com*(seismic_mass+sum(mass_calc(:,2)))- sum(mass_calc(:,4)))/seismic_mass;
    
    mass_tc= 0;
    mass_rc= 0;
    
    mass_node_data= zeros(num_story,6);
    
    mass_node_data(:,1)= [1:num_story]'*10000;
    mass_node_data(:,2)= slab_x_com;
    mass_node_data(:,3)= slab_y_com;
    mass_node_data(:,4)= [story_height:story_height:num_story*story_height]';
    mass_node_data(:,5)= mass_tc;
    mass_node_data(:,6)= mass_rc;
    
    fprintf(key_file, '$Seismic Mass Node\n');
    
    for i =1:size(mass_node_data,1)
        for ii=1:size(mass_node_data,2)
            if ii== size(mass_node_data,2)
                fprintf(key_file, '%d', mass_node_data(i,ii));
            elseif ii< size(mass_node_data,2)
                fprintf(key_file, '%d,', mass_node_data(i,ii));
            end
        end
        fprintf(key_file, '\n');
    end
    
    %*SECTION_BEAM
    %prompt = 'Which element formulation?\nEnter:\n1 for Hughes-Liu\n2 for resultant\n';
    %element_form = input(prompt);
    
    element_form= 2;
    
    section_data1= zeros(size(wall,1),7);
    section_data2= zeros(size(wall,1),6);
    
    section_data1(:,1)=[1:size(wall,1)]';   %section ID
    if element_form==1
        section_data1(:,2)= 1;              %Hughes-Liu beam
    elseif element_form==2
        section_data1(:,2)= 2;              %resultant beam
    else
        printf('\nError\n');
    end
    section_data1(:,3)= 5/6;                %shear factor for rectangle
    section_data1(:,4)= 2;                  %default quandrature rule
    section_data1(:,5)= 0;                  %rectangular section
    section_data1(:,6)= 0;                  %SCOOR
    section_data1(:,7)= 0;                  %non-structural mass per unit length
    
    %s axis: positive to the north
    %t axis: positive to the east
    if element_form==1
        for i=1:size(wall,1)
            if strcmp(wall_orientation(i,1),'E/W')==1
                section_data2(i,1)= wall(i,5);              %beam thickness in s direction at node 1: wall thickness
                section_data2(i,3)= wall(i,4);              %beam thickness in s direction at node 2: wall length
                section_data2(i,2)= section_data2(i,1);     %beam thickness in t direction at node 1: wall thickness
                section_data2(i,4)= section_data2(i,3);     %beam thickness in t direction at node 2: wall length
            elseif strcmp(wall_orientation(i,1),'N/S')==1
                section_data2(i,1)= wall(i,4);              %beam thickness in s direction at node 1: wall length
                section_data2(i,3)= wall(i,5);              %beam thickness in s direction at node 2: wall thickness
                section_data2(i,2)= section_data2(i,1);     %beam thickness in t direction at node 1: wall length
                section_data2(i,4)= section_data2(i,3);     %beam thickness in t direction at node 2: wall thickness
            else
            end
        end
    elseif element_form==2
        %moment of inertia 70% of actual value to account for cracking
        for i=1:size(wall,1)
            section_data2(i,1)= wall(i,4)*wall(i,5);        %area
            section_data2(i,4)= wall(i,4)*wall(i,5)^3/3;    %J
            section_data2(i,5)= section_data2(i,1);         %shear area
            section_data2(i,6)= 0;                          %product moment of inertia

            if strcmp(wall_orientation(i,1),'E/W')==1
                section_data2(i,2)= wall(i,5)*wall(i,4)^3/12*0.7;     %I about s-axis
                section_data2(i,3)= wall(i,4)*wall(i,5)^3/12*0.7;     %I about t-axis
            elseif strcmp(wall_orientation(i,1),'N/S')==1
                section_data2(i,2)= wall(i,4)*wall(i,5)^3/12*0.7;     %I about s-axis
                section_data2(i,3)= wall(i,5)*wall(i,4)^3/12*0.7;     %I about t-axis
            else
            end
        end
    else
    end
    
    fprintf(key_file,'*SECTION_BEAM\n');
    
    for i =1:size(section_data1,1)
        %print section_data1 values
        for ii=1:size(section_data1,2)
            if ii== size(section_data1,2)
                fprintf(key_file, '%.4g', section_data1(i,ii));
                fprintf(key_file, '\n');
            elseif ii< size(section_data1,2)
                fprintf(key_file, '%.4g,', section_data1(i,ii));
            else
            end
        end
        %print section_data2 values
        for ii=1:size(section_data2,2)
            if ii== size(section_data2,2)
                fprintf(key_file, '%.4g', section_data2(i,ii));
                fprintf(key_file, '\n');
            elseif ii< size(section_data2,2)
                fprintf(key_file, '%.4g,', section_data2(i,ii));
            else
            end
        end
    end
        
    %*PART
    part_data= zeros(size(section_data1,1),8);
    
    part_data(:,1)= [1:size(wall,1)]';  %part ID
    part_data(:,2)= [1:size(wall,1)]';  %section ID
    part_data(:,3)= 1;                  %material ID

    fprintf(key_file,'*PART\n$Part ID,Section ID,Material ID\n');

    for i =1:size(part_data,1)
        fprintf(key_file, 'Wall %d\n', part_data(i,1));
        for ii=1:size(part_data,2)
            if ii== size(part_data,2)
                fprintf(key_file, '%d', part_data(i,ii));
            elseif ii< size(part_data,2)
                fprintf(key_file, '%d,', part_data(i,ii));
            end
        end
        fprintf(key_file, '\n');
    end

    %*ELEMENT_BEAM
    element_data= zeros(size(wall,1)*num_story,10);
    
    row= 1;
    
    for story=1:num_story
        element_data(row:row+size(wall,1)-1,1)= story*1000+[1:size(wall,1)]';
        element_data(row:row+size(wall,1)-1,2)= part_data(:,1);
        %node 1 at the top and node 2 at the bottom of each floor
        element_data(row:row+size(wall,1)-1,3)= node_data(row+size(wall,1):row+size(wall,1)*2-1,1);
        element_data(row:row+size(wall,1)-1,4)= node_data(row:row+size(wall,1)-1,1);
        element_data(row:row+size(wall,1)-1,5)= node3_data(:,1);
        
        row= row+ size(wall,1);
    end
    
    %fix all DOF at joints
    element_data(:,6)= 0;
    element_data(:,7)= 0;
    element_data(:,8)= 0;
    element_data(:,9)= 0;
    %default local coordinate system
    element_data(:,10)= 2;

    fprintf(key_file,'*ELEMENT_BEAM\n$Element ID,Part ID,Node 1,Node2,Node3,Node 1 Translational Constraint,Node 2 Translational Constraint,Node 1 Rotational Constraint,Node 2 Rotational Constraint\n');
    
    for i =1:size(element_data,1)
        for ii=1:size(element_data,2)
            if ii== size(element_data,2)
                fprintf(key_file, '%d', element_data(i,ii));
            elseif ii< size(element_data,2)
                fprintf(key_file, '%d,', element_data(i,ii));
            end
        end
        fprintf(key_file, '\n');
    end
    
    %*ELEMENT_MASS
    mass_data= zeros(num_story,3);
    
    mass_data(:,1)= 10000*[1:num_story]';    %element ID
    mass_data(:,2)= mass_data(:,1);         %node of assignment
    mass_data(:,3)= seismic_mass;           %mass value
    
    fprintf(key_file,'*ELEMENT_MASS\n$Element ID,Node of Assignment,Node 1,Mass\n');
    
    for i =1:size(mass_data,1)
        for ii=1:size(mass_data,2)
            if ii== size(mass_data,2)
                fprintf(key_file, '%d', mass_data(i,ii));
            elseif ii< size(mass_data,2)
                fprintf(key_file, '%d,', mass_data(i,ii));
            end
        end
        fprintf(key_file, '\n');
    end
    
    %*BOUNDARY_SPC_NODE
    constraint_data= zeros(size(wall,1),8);
    
    constraint_data(:,1)= wall(:,1)+ 1000;
    constraint_data(:,3:8)= 1;
    
    fprintf(key_file,'*BOUNDARY_SPC_NODE\n$Node ID,Coordinate System ID,DOF-X,DOF-Y,DOF-Z,RDOF-X,RDOF-Y,RDOF-Z\n');
    
    for i =1:size(constraint_data,1)
        for ii=1:size(constraint_data,2)
            if ii== size(constraint_data,2)
                fprintf(key_file, '%d', constraint_data(i,ii));
            elseif ii< size(constraint_data,2)
                fprintf(key_file, '%d,', constraint_data(i,ii));
            end
        end
        fprintf(key_file, '\n');
    end
    
    %*SET_NODE_LIST
    nodeset_data= zeros(num_story,2);
    
    nodeset_data(:,1)= mass_node_data(:,1);
        
    for i=1+size(wall,1):size(node_data,1)
        nodeset_data(floor(node_data(i,1)/1000)-1,mod(node_data(i,1),1000)+1)= node_data(i,1);
    end
        
    for i =1:size(nodeset_data,1)
        node_count= 1;

        for ii=1:size(nodeset_data,2)
            if ii==1
                fprintf(key_file, '*SET_NODE_LIST\n%d\n', i);
            end
            if ii== size(nodeset_data,2)
                fprintf(key_file, '%d', nodeset_data(i,ii));
            elseif ii< size(nodeset_data,2)&& node_count< 8
                fprintf(key_file, '%d,', nodeset_data(i,ii));
                node_count= node_count+ 1;
            %skip to new line when 8 nodes are printed already
            elseif ii< size(nodeset_data,2)&& node_count>= 8
                fprintf(key_file, '%d\n', nodeset_data(i,ii));
                node_count= 1;
            end
        end
        fprintf(key_file, '\n');
    end
    
    %*CONSTRAINED_NODAL_RIGID_BODY
    rigidbody_data= zeros(size(nodeset_data,1),7);
    
    rigidbody_data(:,2)= 0;                     %coorinate system ID
    rigidbody_data(:,4)= 0;                     %output data node
    rigidbody_data(:,5)= 1;                     %print options
    rigidbody_data(:,6)= -3;                    %release z displacement
    rigidbody_data(:,7)= -4;                    %release x and y rotation
    
    for i=1:size(rigidbody_data,1)
        rigidbody_data(i,1)= size(part_data,1)+i;   %part ID
        rigidbody_data(i,3)= i;                     %node set ID
    end
    
    fprintf(key_file,'*CONSTRAINED_NODAL_RIGID_BODY\n');
    
    for i =1:size(rigidbody_data,1)
        for ii=1:size(rigidbody_data,2)
            if ii== size(rigidbody_data,2)
                fprintf(key_file, '%d', rigidbody_data(i,ii));
            elseif ii< size(rigidbody_data,2)
                fprintf(key_file, '%d,', rigidbody_data(i,ii));
            end
        end
        fprintf(key_file, '\n');
    end
    
    %*DEFINE_CURVE
    unit_load_curve= zeros(3,2);
    
    unit_load_curve(1,:)= [0 0];
    unit_load_curve(2,:)= [1 1];
    unit_load_curve(3,:)= [100 1];
    
    loadcurve_data= zeros(size(wall,1),8);
    
    loadcurve_data(:,1)= [1:size(wall,1)]';     %load curve ID
    loadcurve_data(:,2)= 0;                     %dynamic relaxation flag
    loadcurve_data(:,3)= 1;                     %scale factor for x values
    loadcurve_data(:,5)= 0;                     %offset for x values
    loadcurve_data(:,6)= 0;                     %offset for y values
    loadcurve_data(:,7)= 0;                     %data type, 0 for general time dependent curve
    loadcurve_data(:,8)= 0;                     %number of discretization intervals, 100 is dfault
    
    total_area= sum(wall(:,4).*wall(:,5));
    for i=1:size(loadcurve_data,1)
        loadcurve_data(i,4)= wall(i,4)*wall(i,5)/total_area*seismic_mass*9.81;
    end
    
    for i =1:size(loadcurve_data,1)
        fprintf(key_file, '*DEFINE_CURVE\n');
        %print card 1
        for ii=1:size(loadcurve_data,2)
            if ii== size(loadcurve_data,2)
                fprintf(key_file, '%.4g', loadcurve_data(i,ii));
            elseif ii< size(loadcurve_data,2)
                fprintf(key_file, '%.4g,', loadcurve_data(i,ii));
            end
        end
        fprintf(key_file, '\n');
        %print card 2
        for load_curve_row=1:size(unit_load_curve,1)
            for load_curve_column=1:size(unit_load_curve,2)
                if load_curve_column== size(unit_load_curve,2)
                    fprintf(key_file, '%d', unit_load_curve(load_curve_row,load_curve_column));
                elseif load_curve_column< size(unit_load_curve,2)
                    fprintf(key_file, '%d,', unit_load_curve(load_curve_row,load_curve_column));
                end
            end
            fprintf(key_file, '\n');
        end
    end
    
    %*LOAD_NODE
    nodeload_data= zeros(size(wall,1),9);
    
    nodeload_data(:,1)= [1:size(wall,1)]';      %wall ID
    %note that first columnn of  the actual card requres node ID
    nodeload_data(:,3)= 3;                      %direction of load
    nodeload_data(:,4)= [1:size(wall,1)]';      %load curve ID
    nodeload_data(:,5)= -1;                     %load curve scale factor
    nodeload_data(:,6)= 0;                      %coordinate system
    nodeload_data(:,7)= 0;                      %node 1 ID if DOF is a follower force or moment
    nodeload_data(:,8)= 0;                      %node 2 ID if DOF is a follower force or moment
    nodeload_data(:,9)= 0;                      %node 3 ID if DOF is a follower force or moment
    
    fprintf(key_file, '*LOAD_NODE_POINT\n');
    
    for story_index=1:num_story
        nodeload_data(:,2)= nodeload_data(:,1)+1000*(story_index+1);
        for i =1:size(nodeload_data,1)
            for ii=2:size(nodeload_data,2)
                if ii== size(nodeload_data,2)
                    fprintf(key_file, '%d', nodeload_data(i,ii));
                elseif ii< size(nodeload_data,2)
                    fprintf(key_file, '%d,', nodeload_data(i,ii));
                end
            end
            fprintf(key_file, '\n');
        end
    end
    
    %*CONTROL_IMPLICIT_SOLVER
    fprintf(key_file, '*CONTROL_IMPLICIT_SOLVER\n');
    fprintf(key_file, ',3\n');
    
    %Eigenvalue problem
    fprintf(key_file, '$ Note: Need to adjust boundary conditions to retrain movement out of plane.\n');
    fprintf(key_file, '*CONTROL_IMPLICIT_GENERAL\n');
    fprintf(key_file, '$IMFLAG	DT0\n');
    fprintf(key_file, '1,5.00E-03\n');
    fprintf(key_file, '*CONTROL_IMPLICIT_EIGENVALUE\n');
    fprintf(key_file, '$NEIG	CENTER	IFLAG	LFTEND	RFLAG	RHTEND	EIGMTH	SHFSCL\n');
    fprintf(key_file, '10,1\n');
    fprintf(key_file, '$ISOLID	IBEAM	ISHELL	ITSHELL	MSTRESS	EVDUMP\n');
    fprintf(key_file, '0,0,0,0,0,-1\n');
    
    %*END
    fprintf(key_file,'*END');
    fclose(key_file);
    
    %create folder
    folder_name= sprintf('bldg%d_file',bldg_id);
    mkdir(folder_name)
    
    %move file to new folder
    movefile(sprintf('Bldg%d_Keyword_File.key',bldg_id),folder_name)

end

