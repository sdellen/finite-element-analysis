clear; clc;

% Q1: ----------------
% joints = [0, 0;
%          0.08, 0;
%          0.24, 0;
%          0.48, 0];
% 
% members = [1, 2, 69e9, (pi*1.5^2)*1e-4;
%             2, 3, 69e9, (pi*1^2)*1e-4;
%             3, 4, 69e9, (pi*.5^2)*1e-4];
% 
% zero_dx = [1];
% zero_dy = [1, 2, 3, 4];
% 
% applied_fx = [4, 6100];
% applied_fy = [];
% 
% joint_displacements = do_FEA(joints, members, zero_dx, zero_dy, applied_fx, applied_fy);
% Q1: ----------------

% Q2: ----------------

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

joint_displacements = do_FEA(joints, members, zero_dx, zero_dy, applied_fx, applied_fy);
% Q2: ----------------

function result = do_FEA(joints, members, zero_dx, zero_dy, applied_fx, applied_fy)

    num_members = size(members,1);
    num_joints = size(joints,1);

    E = members(:, 3);
    A = members(:, 4);
    
    dx = joints(members(:, 2), 1) - joints(members(:, 1), 1);
    dy = joints(members(:, 2), 2) - joints(members(:, 1), 2);
    
    L = sqrt(dx.^2 + dy.^2);
    t = atan(dy ./ dx);

    k = E.*A./L;
    
    K = truss_stiffness_g(k, t, members, num_members, num_joints);
    
    dl = K;
    
    relevent_forces = zeros(num_joints*2, 1);
    
    for i = 1:size(applied_fx, 1)
        relevent_forces(max(applied_fx(i, 1)*2-1, 1), 1) = applied_fx(i, 2);
    end
    
    for i = 1:size(applied_fy, 1)
        relevent_forces(max(applied_fy(i, 1)*2, 2), 1) = applied_fy(i, 2);
    end
    
    relevent_forces = relevent_forces(~ismember(1:num_joints * 2, union(zero_dx * 2 - 1, zero_dy * 2)));
    
    for i = 1:length(zero_dx)
        dl(:,max(zero_dx(i)*2-1, 1)) = 0;
        dl(max(zero_dx(i)*2-1, 1),:) = 0;
    end
    
    for i = 1:length(zero_dy)
        dl(:,max(zero_dy(i)*2, 2)) = 0;
        dl(max(zero_dy(i)*2, 2),:) = 0;
    end
    
    dl = dl(any(dl), any(dl));
    dl = dl \ relevent_forces;

    dl_temp = dl;
    d = zeros(num_joints * 2, 1);

    for i = 1:num_joints*2
        if mod(i, 2) == 1
            if ~ismember((i+1) / 2, zero_dx)
                d(i) = dl_temp(1);
                dl_temp = dl_temp(2:end);
            end
        else
            if ~ismember(i / 2, zero_dy)
                d(i) = dl_temp(1);
                dl_temp = dl_temp(2:end);
            end
        end
    end
    
    f = K * d;
    
    s = zeros(length(E), 1);
    
    for i = 1:length(s)
        fn = 2*members(i, 1)-1;
        sn = 2*members(i, 2)-1;
    
        s(i, 1) = E(i) / L(i) * [-1, 1] * [cos(t(i)) sin(t(i)) 0 0;0 0 cos(t(i)) sin(t(i))] * [d(fn:fn+1);d(sn:sn+1)];
    end

    fprintf("------- Total Stiffness Matrix\n\n");
    disp(K);
    
    fprintf("------- Joint Forces:\n\n");
    for i = 1:length(f)/2
        fprintf(format_number("F"+i+"x", f(2*i-1), 'N'))
        fprintf(format_number("F"+i+"y", f(2*i), 'N'))

    end
    
    fprintf("\n------- Joint Displacements:\n\n");
    for i = 1:length(f)/2
        fprintf(format_number("d"+i+"x", d(2*i-1), 'm'))
        fprintf(format_number("d"+i+"y", d(2*i), 'm'))
    end
    
    fprintf("\n------- Member Stresses:\n\n");
    for i = 1:length(s)
        fprintf(format_number("S"+i, s(i), 'Pa'))
    end
    
    
    fprintf("\n-------\n");

    result = d;

end

function result = truss_stiffness_g(k, t, members, num_members, num_joints)
    if length(k) ~= length(t)
        error('The two matrices must be of the same length');
    end

    K = zeros(2*num_joints, 2*num_joints);

    for i = 1:num_members
        C = cos(t(i, 1));
        S = sin(t(i, 1));
        kb = k(i) * [C^2, C*S,-C^2,-C*S;
                     C*S, S^2,-C*S,-S^2;
                    -C^2,-C*S, C^2, C*S;
                    -C*S,-S^2, C*S, S^2];

        fn = 2*members(i, 1) - 1;
        sn = 2*members(i, 2) - 1;

        K(fn:fn+1,fn:fn+1) = K(fn:fn+1,fn:fn+1) + kb(1:2,1:2);
        K(fn:fn+1,sn:sn+1) = K(fn:fn+1,sn:sn+1) + kb(1:2,3:4);
        K(sn:sn+1,fn:fn+1) = K(sn:sn+1,fn:fn+1) + kb(3:4,1:2);
        K(sn:sn+1,sn:sn+1) = K(sn:sn+1,sn:sn+1) + kb(3:4,3:4);
    end
    result = K;
end


function formatted_string = format_number(variable_name, number, unit)

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