clc; clear; close all;

% Duomenys
n = 5;                          % Virsuniu skaicius
V = 1:n;                        % Virsuniu sarasas
U = [1 2; 2 3; 3 1; 4 1; 1 5]; % Briaunu matrica
m = size(U, 1);                 % Briaunu skaicius

% Jungumo komponentės - MANUAL IMPLEMENTATION
JungiosiosKomp = zeros(1, n);
JungiujuKSk = 0;

while any(JungiosiosKomp == 0)
    NeturinciosJungK = V(JungiosiosKomp == 0);
    PirmaBeJungK = NeturinciosJungK(1);

    JungiujuKSk = JungiujuKSk + 1;
    JungiosiosKomp(PirmaBeJungK) = JungiujuKSk;

    ArTesti = true;
    while ArTesti
        JungiosiosKompBackup = JungiosiosKomp;
        for BriaunosNr = 1:m
            if (JungiosiosKomp(U(BriaunosNr,1)) == JungiujuKSk)
                JungiosiosKomp(U(BriaunosNr,2)) = JungiujuKSk;
            end
            if (JungiosiosKomp(U(BriaunosNr,2)) == JungiujuKSk)
                JungiosiosKomp(U(BriaunosNr,1)) = JungiujuKSk;
            end
        end
        ArTesti = any(JungiosiosKompBackup ~= JungiosiosKomp);
    end
end

disp('JungiosiosKomp:');
disp(JungiosiosKomp);

% Randame medzius
BriaunuKiek = n - 1;
RastiMedziai = []; 
RastiMedziaiCount = 0;

% Manual combinations generator
function combos = generateCombinations(arr, k)
    if k == 0
        combos = zeros(1,0);
        return
    elseif k == length(arr)
        combos = arr;
        return
    end
    
    if k == 1
        combos = arr';
        return
    end
    
    combos = [];
    for i = 1:length(arr)-k+1
        subcombos = generateCombinations(arr(i+1:end), k-1);
        for j = 1:size(subcombos,1)
            combos = [combos; arr(i), subcombos(j,:)]; %#ok<AGROW>
        end
    end
end

if m >= BriaunuKiek
    GalimiMedziai = generateCombinations(1:m, BriaunuKiek); % Manual combinations
    for i = 1:size(GalimiMedziai,1)
        Eil = GalimiMedziai(i,:);
        Briaunos = U(Eil,:);

        % Create adjacency matrix
        Adj = zeros(n,n);
        for j = 1:size(Briaunos,1)
            u = Briaunos(j,1);
            v = Briaunos(j,2);
            Adj(u,v) = 1;
            Adj(v,u) = 1;
        end

        % Check connectivity manually
        visited = zeros(1,n);
        queue = 1;
        visited(1) = 1;

        while ~isempty(queue)
            current = queue(1);
            queue(1) = [];
            neighbors = find(Adj(current,:) > 0);
            for neighbor = neighbors
                if visited(neighbor) == 0
                    visited(neighbor) = 1;
                    queue(end+1) = neighbor;
                end
            end
        end

        if all(visited == 1)
            RastiMedziaiCount = RastiMedziaiCount + 1;
            RastiMedziai(RastiMedziaiCount,:,:) = Briaunos; %#ok<SAGROW>
        end
    end
end

% Plot original graph
Ucell = U; % Simple array works
title('Pradinis grafas');
plotGraphVU1(V, num2cell(Ucell,2), 0, 0, [], 0, 10, 1, 'b');

% Plot all spanning trees
for k = 1:size(RastiMedziai,1)
    Medis = squeeze(RastiMedziai(k,:,:));
    figure;
    plotGraphVU1(V, num2cell(Medis,2), 0, 0, [], 0, 10, 1, 'b');
end

disp('Rastų medžių skaičius:');
disp(RastiMedziaiCount);
