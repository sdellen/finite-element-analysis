% clear; clc;

nodes = [0, 0; 4, 0; 0, -2; 2, -2];
elements = [1, 2, 200e9, 0.0004; 2, 4, 200e9, 0.0004; 1, 4, 200e9, 0.0004; 3, 4, 200e9, 0.0004];

known_dx = [1, 3];
known_dy = [1, 3];

applied_fx = [];
applied_fy = [2, -51000];

pins = [1, 3];

num_elements = length(elements);
num_nodes = length(nodes);



E = [elements(1:end, 3)];
A = [elements(1:end, 4)];
L = sqrt((nodes(elements(1:num_elements, 2), 1) -  (nodes(elements(1:num_elements, 1), 1))).^2 + (nodes(elements(1:num_elements, 2), 2) -  (nodes(elements(1:num_elements, 1), 2))).^2);
t = atan((nodes(elements(1:num_elements, 2), 2) -  (nodes(elements(1:num_elements, 1), 2)))./(nodes(elements(1:num_elements, 2), 1) -  (nodes(elements(1:num_elements, 1), 1)))); % angles for elements

k = E.*A./L;

for e = 1:num_elements
    % Extract element properties
    node1 = elements(e, 1);
    node2 = elements(e, 2);
    E(e,1) = elements(e, 3); % Young's modulus
    
    % Nodal coordinates
    x1 = nodes(node1, 1);
    y1 = nodes(node1, 2);
    x2 = nodes(node2, 1);
    y2 = nodes(node2, 2);
    
    % Element length and cosine and sine values
    % L(e,1) = sqrt((x2 - x1)^2 + (y2 - y1)^2);
    % t(1,e) = atan((y2 - y1)/(x2 - x1));
    % c = (x2 - x1) / L;
    % s = (y2 - y1) / L;
end

ke = truss_stiffness(k, t);
K = truss_stiffness_g(k, t, elements, num_elements, num_nodes);

dl = K;

relevent_forces = zeros(num_nodes*2, 1);

% disp(size(applied_fy, 1));
for i = 1:size(applied_fx, 1)
    relevent_forces(max(applied_fx(i, 1)*2-1, 1), 1) = applied_fx(i, 2);
end

for i = 1:size(applied_fy, 1)
    relevent_forces(max(applied_fy(i, 1)*2, 2), 1) = applied_fy(i, 2);
end

for i = 1:length(known_dx)
    dl(:,max(known_dx(i)*2-1, 1)) = zeros(num_nodes*2, 1);
    dl(max(known_dx(i)*2-1, 1),:) = zeros(1, num_nodes*2);
    % relevent_forces(max(known_dx(i)*2-1, 1), :) = [];
end

for i = 1:length(known_dy)
    dl(:,max(known_dy(i)*2, 2)) = zeros(num_nodes*2, 1);
    dl(max(known_dy(i)*2, 2),:) = zeros(1, num_nodes*2);
    % relevent_forces(max(known_dy(i)*2, 2), :) = [];
end

dl = dl(any(dl), any(dl));

disp(dl);
disp(relevent_forces);





dl = [K(3:4,3:4) K(3:4,7:8);K(7:8,3:4) K(7:8,7:8)]\...
    [0;-51000;0;0];

d = [0;0;dl(1:2);0;0;dl(3:end)];

f = K * d;

s = zeros(length(E), 1);

for i = 1:length(s)
    fn = 2*elements(i, 1)-1;
    sn = 2*elements(i, 2)-1;

    s(i, 1) = E(i) / L(i) * [-1, 1] * [cos(t(i)) sin(t(i)) 0 0;0 0 cos(t(i)) sin(t(i))] * [d(fn:fn+1);d(sn:sn+1)];
end

fprintf("Total Stiffness Matrix\n\n");
disp(K);

fprintf("\n-------\nNodal Forces:\n\n");
for i = 1:length(f)/2
    fprintf("F%ix = %7.3f\n", i, f(2*i-1));
    fprintf("F%iy = %7.3f\n", i, f(2*i));
end

fprintf("\n-------\nNodal Displacements:\n\n");
for i = 1:length(f)/2
    fprintf("d%ix = %7.3f\n", i, d(2*i-1));
    fprintf("d%iy = %7.3f\n", i, d(2*i));
end

fprintf("\n-------\nElement Stresses:\n\n");
for i = 1:length(s)
    fprintf("S%i = %8.3f\n", i, s(i));
end


fprintf("\n-------\n");



function result = truss_stiffness(k, t)
    if length(k) ~= length(t)
        error('The two matrices must be of the same length');
    end


    K = zeros(4, 4, length(k));


    for i = 1:length(k)
        C = cos(t(i, 1));
        S = sin(t(i, 1));
        kb = k(i) * [C^2, C*S,-C^2,-C*S;
                     C*S, S^2,-C*S,-S^2;
                    -C^2,-C*S, C^2, C*S;
                    -C*S,-S^2, C*S, S^2];
        
        K(:,:,i) = kb;
    end

    result = K;
    
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
        % disp(K);
    end
    result = K;
end
