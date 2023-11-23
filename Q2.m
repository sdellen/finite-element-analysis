clear; clc;

ult_tensile = 500e6;
ult_compressive = -250e6;

joints = [0, 0;
         4, 0;
         0, -2;
         2, -2];

members = [1, 2, 200e9, 0.0004;
            2, 4, 200e9, 0.0004;
            1, 4, 200e9, 0.0004;
            3, 4, 200e9, 0.0004];

zero_dx = [1, 3];
zero_dy = [1, 3];

applied_fx = [];
applied_fy = [2, -51000];

% Perform FEA on the truss structure
joint_displacements = do_FEA(joints, members, zero_dx, zero_dy, applied_fx, applied_fy);

function result = do_FEA(joints, members, zero_dx, zero_dy, applied_fx, applied_fy)
    %{
    Function to perform Finite Element Analysis (FEA) on a truss structure

    Inputs:
      joints: Matrix containing the coordinates of joints
      members: Matrix containing the connectivity, Young's Modulus, and cross-sectional areas of members
      zero_dx: Vector specifying joints with zero displacement in the x-direction
      zero_dy: Vector specifying joints with zero displacement in the y-direction
      applied_fx: Matrix specifying applied forces at specified joints in the x-direction
      applied_fy: Matrix specifying applied forces at specified joints in the y-direction

    Output:
      result: Nodal displacements in the truss structure after FEA
    %}

    % Extract the number of members and joints
    num_members = size(members, 1);
    num_joints = size(joints, 1);

    % Extract Young's Modulus and cross-sectional area of the members
    E = members(:, 3);
    A = members(:, 4);

    % Calculate the differences in x and y coordinates between joints per member
    dx = joints(members(:, 2), 1) - joints(members(:, 1), 1);
    dy = joints(members(:, 2), 2) - joints(members(:, 1), 2);

    % Calculate the length, angle, and stiffness of each member
    member_lengths = sqrt(dx.^2 + dy.^2);
    member_angles = atan(dy ./ dx);
    k = E .* A ./ member_lengths;

    % Calculate the stiffness matrix for the truss structure
    K = truss_stiffness_g(k, member_angles, members, num_members, num_joints);

    dl = K;
    relevent_forces = zeros(num_joints * 2, 1);

    % Apply external forces at specified joints
    for i = 1:size(applied_fx, 1)
        relevent_forces(max(applied_fx(i, 1) * 2 - 1, 1), 1) = applied_fx(i, 2);
    end

    for i = 1:size(applied_fy, 1)
        relevent_forces(max(applied_fy(i, 1) * 2, 2), 1) = applied_fy(i, 2);
    end

    % Remove constrained degrees of freedom from the force vector
    relevent_forces = relevent_forces(~ismember(1:num_joints * 2, union(zero_dx * 2 - 1, zero_dy * 2)));

    % set displacements known to be constrained to 0
    for i = 1:length(zero_dx)
        dl(:, max(zero_dx(i) * 2 - 1, 1)) = 0;
        dl(max(zero_dx(i) * 2 - 1, 1), :) = 0;
    end

    for i = 1:length(zero_dy)
        dl(:, max(zero_dy(i) * 2, 2)) = 0;
        dl(max(zero_dy(i) * 2, 2), :) = 0;
    end

    % Solve for the nodal displacements using the stiffness matrix
    dl = dl(any(dl), any(dl));
    dl = dl \ relevent_forces;

    % Reconstruct the full displacement vector with zero displacements
    d = zeros(num_joints * 2, 1);

    % Populate the displacement vector with non-zero displacements
    for i = 1:num_joints * 2
        if mod(i, 2) == 1
            if ~ismember((i + 1) / 2, zero_dx)
                d(i) = dl(1);
                dl = dl(2:end);
            end
        else
            if ~ismember(i / 2, zero_dy)
                d(i) = dl(1);
                dl = dl(2:end);
            end
        end
    end

    % Calculate nodal forces in each member
    f = K * d;

    % Calculate stresses in each member
    s = zeros(length(E), 1);

    for i = 1:length(s)
        fn = 2 * members(i, 1) - 1;
        sn = 2 * members(i, 2) - 1;

        s(i, 1) = E(i) / member_lengths(i) * [-1, 1] * [cos(member_angles(i)), sin(member_angles(i)), 0, 0; 0, 0, cos(member_angles(i)), sin(member_angles(i))] * [d(fn:fn+1); d(sn:sn+1)];
    end

    % Display results
    fprintf("------- Total Stiffness Matrix\n\n");
    disp(K);

    fprintf("------- Joint Forces:\n\n");
    for i = 1:length(f) / 2
        fprintf(format_number("F" + i + "x", f(2 * i - 1), 'N'))
        fprintf(format_number("F" + i + "y", f(2 * i), 'N'))

    end

    fprintf("\n------- Joint Displacements:\n\n");
    for i = 1:length(f) / 2
        fprintf(format_number("d" + i + "x", d(2 * i - 1), 'm'))
        fprintf(format_number("d" + i + "y", d(2 * i), 'm'))
    end

    fprintf("\n------- Member Stresses:\n\n");
    for i = 1:length(s)
        fprintf(format_number("S" + i, s(i), 'Pa'))
    end

    fprintf("\n-------\n");

    % Return the nodal displacements as the result
    result = d;
end

function result = truss_stiffness_g(k, member_angle, members, num_members, num_joints)
    %{
    Function to calculate the global stiffness matrix for a truss structure
    
    Inputs:
      k: Vector of member stiffness values
      t: Vector of member angles (in radians)
      members: Matrix defining the connectivity of members
      num_members: Number of members in the truss
      num_joints: Number of joints in the truss

    Output:
      result: Global stiffness matrix for the entire truss
    %}
    
    % Check if the input vectors have the same length
    if length(k) ~= length(member_angle)
        error('The two matrices must be of the same length');
    end

    % Initialize the global stiffness matrix
    K = zeros(2 * num_joints, 2 * num_joints);

    % Loop through each member to assemble the stiffness matrix
    for i = 1:num_members
        % Calculate cosine and sine of the member angle
        C = cos(member_angle(i, 1));
        S = sin(member_angle(i, 1));

        % Calculate the local stiffness matrix for the current member
        kb = k(i) * [C^2, C*S, -C^2, -C*S;
                     C*S, S^2, -C*S, -S^2;
                    -C^2, -C*S, C^2, C*S;
                    -C*S, -S^2, C*S, S^2];

        % Define the degrees of freedom for the start and end joints of the member
        fn = 2 * members(i, 1) - 1;
        sn = 2 * members(i, 2) - 1;

        % Assemble the local stiffness matrix into the global stiffness matrix
        K(fn:fn+1, fn:fn+1) = K(fn:fn+1, fn:fn+1) + kb(1:2, 1:2);
        K(fn:fn+1, sn:sn+1) = K(fn:fn+1, sn:sn+1) + kb(1:2, 3:4);
        K(sn:sn+1, fn:fn+1) = K(sn:sn+1, fn:fn+1) + kb(3:4, 1:2);
        K(sn:sn+1, sn:sn+1) = K(sn:sn+1, sn:sn+1) + kb(3:4, 3:4);
    end

    % Output the global stiffness matrix as the result
    result = K;
end

function formatted_string = format_number(variable_name, number, unit)
    %{
    Function to format numbers with appropriate prefixes and units

    Inputs:
      variable_name: Name of the variable or quantity
      number: Numeric value to be formatted
      unit: Unit of the quantity

    Output:
      formatted_string: Formatted string representing the variable, number, and unit
    %}
    
    % Check if the number is zero
    if number == 0
        % No suffix for 0
        formatted_string = sprintf('%s = %10.5f %s\n', variable_name, number, unit);
    else
        % Extract the sign of the number
        sign = true;
        if number < 0
            sign = false;
            number = abs(number);
        end
    
        % Define small and big suffixes
        small_suffixes = {'', 'm', 'μ', 'n', 'p', 'f', 'a', 'z', 'y'};
        big_suffixes = {'', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'};
    
        % Find the appropriate suffix and scale factor
        exp_val = floor(log10(number));
        prefix_id = min(floor(abs(exp_val) / 3) + 1 + ((number < 1) * (1 - (mod(exp_val, 3) == 0))), 9);
        scale_factor = 10^( (1-2*(number < 1)) * 3 * (prefix_id - 1) );
    
        % Format the number
        formatted_number = sprintf('%10.5f', (2*sign-1) * number / scale_factor);

        % Construct the formatted string
        formatted_string = sprintf('%s = %s %s%s\n', variable_name, formatted_number, (number < 1) * small_suffixes{prefix_id} + (number >= 1) * big_suffixes{prefix_id}, unit);
    end
end
