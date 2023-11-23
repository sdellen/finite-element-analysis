clear; clc;


% Q1: ----------------
nodes = [0, 0;
         0.08, 0;
         0.24, 0;
         0.48, 0];


elements = [1, 2, 69e9, (pi*1.5^2)*1e-4;
            2, 3, 69e9, (pi*1^2)*1e-4;
            3, 4, 69e9, (pi*.5^2)*1e-4];

zero_dx = [1];
zero_dy = [1, 2, 3, 4];

applied_fx = [4, 6100];
applied_fy = [];

nodal_displacements = do_FEA(nodes, elements, zero_dx, zero_dy, applied_fx, applied_fy);
% Q1: ----------------

% Q2: ----------------
% nodes = [0, 0; 4, 0; 0, -2; 2, -2];
% elements = [1, 2, 200e9, 0.0004; 2, 4, 200e9, 0.0004; 1, 4, 200e9, 0.0004; 3, 4, 200e9, 0.0004];
% 
% zero_dx = [1, 3];
% zero_dy = [1, 3];
% 
% applied_fx = [];
% applied_fy = [2, -51000];
% 
% nodal_displacements = do_FEA(nodes, elements, zero_dx, zero_dy, applied_fx, applied_fy);
% Q2: ----------------

function result = do_FEA(nodes, elements, zero_dx, zero_dy, applied_fx, applied_fy)

    num_elements = size(elements,1);
    num_nodes = size(nodes,1);
    E = [elements(1:end, 3)];
    A = [elements(1:end, 4)];
    L = sqrt((nodes(elements(1:num_elements, 2), 1) -  (nodes(elements(1:num_elements, 1), 1))).^2 + (nodes(elements(1:num_elements, 2), 2) -  (nodes(elements(1:num_elements, 1), 2))).^2);
    t = atan((nodes(elements(1:num_elements, 2), 2) -  (nodes(elements(1:num_elements, 1), 2)))./(nodes(elements(1:num_elements, 2), 1) -  (nodes(elements(1:num_elements, 1), 1)))); % angles for elements
    
    k = E.*A./L;
    
    K = truss_stiffness_g(k, t, elements, num_elements, num_nodes);
    
    dl = K;
    
    relevent_forces_temp = zeros(num_nodes*2, 1);
    relevent_forces = [];
    
    for i = 1:size(applied_fx, 1)
        relevent_forces_temp(max(applied_fx(i, 1)*2-1, 1), 1) = applied_fx(i, 2);
    end
    
    for i = 1:size(applied_fy, 1)
        relevent_forces_temp(max(applied_fy(i, 1)*2, 2), 1) = applied_fy(i, 2);
    end
    
    for i = 1:size(relevent_forces_temp)
        if mod(i, 2) == 0
            if zero_dy(:) ~= i / 2
                relevent_forces = [relevent_forces; relevent_forces_temp(i)];
            end
        else
            if zero_dx(:) ~= (i+1) / 2
                relevent_forces = [relevent_forces; relevent_forces_temp(i)];
            end
        end
    end
    
    for i = 1:length(zero_dx)
        dl(:,max(zero_dx(i)*2-1, 1)) = zeros(num_nodes*2, 1);
        dl(max(zero_dx(i)*2-1, 1),:) = zeros(1, num_nodes*2);
    end
    
    for i = 1:length(zero_dy)
        dl(:,max(zero_dy(i)*2, 2)) = zeros(num_nodes*2, 1);
        dl(max(zero_dy(i)*2, 2),:) = zeros(1, num_nodes*2);
    end
    
    dl = dl(any(dl), any(dl));
    dl = dl\relevent_forces;
    dl_temp = dl;
    d = [];

    for i = 1:num_nodes*2
        if mod(i, 2) == 1
            if ismember((i+1) / 2, zero_dx)
                d = [d; 0];
            else
                d = [d; dl_temp(1, 1)];
                dl_temp = dl_temp(2:end, 1);
            end
        else
            if ismember(i / 2, zero_dy)
                d = [d; 0];
            else
                d = [d; dl_temp(1, 1)];
                dl_temp = dl_temp(2:end, 1);
            end
        end
    end
    
    f = K * d;
    
    s = zeros(length(E), 1);
    
    for i = 1:length(s)
        fn = 2*elements(i, 1)-1;
        sn = 2*elements(i, 2)-1;
    
        s(i, 1) = E(i) / L(i) * [-1, 1] * [cos(t(i)) sin(t(i)) 0 0;0 0 cos(t(i)) sin(t(i))] * [d(fn:fn+1);d(sn:sn+1)];
    end
    d = d*1e6;
    fprintf("Total Stiffness Matrix\n\n");
    disp(K);
    
    fprintf("\n-------\nNodal Forces:\n\n");
    for i = 1:length(f)/2
        fprintf("F%ix = %7.5f\n", i, f(2*i-1));
        fprintf("F%iy = %7.5f\n", i, f(2*i));
    end
    
    fprintf("\n-------\nNodal Displacements:\n\n");
    for i = 1:length(f)/2
        fprintf("d%ix = %7.5f\n", i, d(2*i-1));
        fprintf("d%iy = %7.5f\n", i, d(2*i));
    end
    
    fprintf("\n-------\nElement Stresses:\n\n");
    for i = 1:length(s)
        fprintf("S%i = %8.5f\n", i, s(i));
    end
    
    
    fprintf("\n-------\n");

    result = d;

end

function result = truss_stiffness_g(k, t, elements, num_elements, num_nodes)
    if length(k) ~= length(t)
        error('The two matrices must be of the same length');
    end

    K = zeros(2*num_nodes, 2*num_nodes);

    for i = 1:num_elements
        C = cos(t(i, 1));
        S = sin(t(i, 1));
        kb = k(i) * [C^2, C*S,-C^2,-C*S;
                     C*S, S^2,-C*S,-S^2;
                    -C^2,-C*S, C^2, C*S;
                    -C*S,-S^2, C*S, S^2];

        fn = 2*elements(i, 1) - 1;
        sn = 2*elements(i, 2) - 1;

        K(fn:fn+1,fn:fn+1) = K(fn:fn+1,fn:fn+1) + kb(1:2,1:2);
        K(fn:fn+1,sn:sn+1) = K(fn:fn+1,sn:sn+1) + kb(1:2,3:4);
        K(sn:sn+1,fn:fn+1) = K(sn:sn+1,fn:fn+1) + kb(3:4,1:2);
        K(sn:sn+1,sn:sn+1) = K(sn:sn+1,sn:sn+1) + kb(3:4,3:4);
    end
    result = K;
end
