clc; clear; close all;

% Duomenys
n = 5;                          % Virsuniu skaicius
V = 1:n;                        % Virsuniu sarasas
U = [1 2; 2 3; 3 1; 4 1; 1 5];       % Briaunu matrica
m = length(U(:, 1));            % Briaunu skaicius

Ucell = [];
for BriaunosNr = m:-1:1
    Ucell{BriaunosNr} = U(BriaunosNr, :);
end

% Trukmes matavimo pradzia
tic;

% Jungumo komponentės - MANUAL IMPLEMENTATION (without conncomp)
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
            if (JungiosiosKomp(U(BriaunosNr, 1)) == JungiujuKSk)
                JungiosiosKomp(U(BriaunosNr, 2)) = JungiujuKSk;
            end
            if (JungiosiosKomp(U(BriaunosNr, 2)) == JungiujuKSk)
                JungiosiosKomp(U(BriaunosNr, 1)) = JungiujuKSk;
            end
        end
        ArTesti = any(JungiosiosKompBackup ~= JungiosiosKomp);
    end
end

disp('JungiosiosKomp:');
disp(JungiosiosKomp);

% Randame medzius
BriaunuKiek = n - 1;
RastiMedziai = {}; 

if m >= BriaunuKiek
    GalimiMedziai = nchoosek(1:m, BriaunuKiek); % Randa skaiciu kombinacijas (medzius) is skaiciu nuo 1 iki m, elementu kiekis - briaunuKiek
    for i = 1:size(GalimiMedziai,1) % Tikrina kiekviena eilute GalimiMedziai
        Eil = GalimiMedziai(i,:);   % Paima i eilute
        Briaunos = U(Eil,:);        % Paima briaunas is U kuriu numeriai yra isvardinti eil
        
        % CHECK IF GRAPH IS CONNECTED - MANUAL IMPLEMENTATION
        % Create adjacency matrix for this tree candidate
        Adj = zeros(n, n);
        for j = 1:size(Briaunos, 1)
            u = Briaunos(j, 1);
            v = Briaunos(j, 2);
            Adj(u, v) = 1;
            Adj(v, u) = 1;
        end
        
        % Check connectivity using manual BFS
        visited = zeros(1, n);
        queue = 1; % Start from vertex 1
        visited(1) = 1;
        
        while ~isempty(queue)
            current = queue(1);
            queue(1) = [];
            
            % Find all neighbors
            neighbors = find(Adj(current, :) > 0);
            for neighbor = neighbors
                if visited(neighbor) == 0
                    visited(neighbor) = 1;
                    queue(end+1) = neighbor;
                end
            end
        end
        
        % If all vertices are visited, the graph is connected
        if all(visited == 1)
            RastiMedziai{end+1} = Briaunos;   % Jeigu grafas jungusis, pridedame prie RastiMedziai
        end
    end
end

title('Pradinis grafas');
plotGraphVU1(V, Ucell, 0, 0, [], 0, 10, 1, 'b');

for k = 1:length(RastiMedziai) % Einame pro visus rastus medzius
    MedisK = RastiMedziai{k}; % Paima k medi
    Medis = cell(1, size(MedisK,1));
    for j = 1:size(MedisK,1)
        Medis{j} = MedisK(j, :);
    end
    figure;
    plotGraphVU1(V, Medis, 0, 0, [], 0, 10, 1, 'b');
end

% Trukmes matavimo pabaiga
Trukme = toc;
disp('Skaiciavimu trukme, sekundemis:');
disp(Trukme);

disp('Rastų medžių skaičius:');
disp(numel(RastiMedziai));