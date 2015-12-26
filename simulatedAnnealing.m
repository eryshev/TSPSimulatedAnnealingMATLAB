function [result] = simulatedAnnealing(cities, initialTemperature, endTemperature)
n = size(cities,1); % �������� ������ ������� �������

state = randperm(n); % ����� ��������� ���������, ��� ��������� �������
% ������� randperm(n) - ���������� ���������� ������������������ �� ����� ����� �� 1 �� n
delete(findall(0,'Type','figure'))
scatter(cities(:, 1), cities(:, 2));
hold on;
for i = 1:n-1
    line([cities(state(i), 1), cities(state(i+1), 1)], [cities(state(i), 2), cities(state(i+1), 2)]);
end
line([cities(state(1), 1), cities(state(end), 1)], [cities(state(1), 2), cities(state(end), 2)]);
hold off;

currentEnergy = calculateEnergy(state, cities); % ��������� ������� ��� ������� ���������
currentEnergy
T = initialTemperature;

for i = 1:10000  %�� ������ ������ ������������ ���������� ��������
    % ����� ���� ������� ��� ������������ ������� ������� ��������� ����������� T 
    stateCandidate = generateStateCandidate(state.').'; % �������� ���������-��������
    candidateEnergy = calculateEnergy(stateCandidate, cities); % ��������� ��� �������
    
    if(candidateEnergy < currentEnergy) % ���� �������� �������� ������� ��������
        currentEnergy = candidateEnergy; % �� ��� ���������� ������� ����������
        state = stateCandidate;
    else
        p = getTransitionProbability(candidateEnergy-currentEnergy, T); % �����, ������� �����������
        if (isTransition(p)) % � �������, ������������ �� �������
            currentEnergy = candidateEnergy;
            state = stateCandidate;
        end
    end;
    
    T = decreaseTemperature(initialTemperature, i); % ��������� �����������
    if(T <= endTemperature) % ������� ������
        break;
    end;
end

currentEnergy
figure
scatter(cities(:, 1), cities(:, 2));
hold on;
for i = 1:n-1
    line([cities(state(i), 1), cities(state(i+1), 1)], [cities(state(i), 2), cities(state(i+1), 2)]);
end
line([cities(state(1), 1), cities(state(end), 1)], [cities(state(1), 2), cities(state(end), 2)]);
hold off;
result = state;
end


function [ T ] = decreaseTemperature( initialTemperature, i)
%initialTemperature - ��������� �����������
%i - ����� ��������
T = initialTemperature * 0.1 / i;
end

function [ P ] = getTransitionProbability( dE, T )
P = exp(-dE/T);
end

function [ a ] = isTransition(probability)
value = rand(1);

if(value <= probability)
    a = 1;
else
    a = 0;
end
end

function [ res ] = generateStateCandidate(seq)
%seq - ���������� ��������� (�������), �� ��������
%�� ����� �������� ���������-��������

n = size(seq,1); % ���������� ������ ������������������
l = randi(n,1); % ���������� ����� ��������� �����
r = randi(n,1);% ���������� ����� ��������� �����
res = seq;
if(l > r)
    res(r:l) = flipud(seq(r:l)); % �������� ���������������������
else
    res(l:r) = flipud(seq(l:r)); % �������� ���������������������
end

end

function [ E ] = calculateEnergy(sequence, cities)
    function [ distance ] = metric( A, B )
        distance = (A - B).^2;
        distance = sqrt(distance);
        distance = sum(distance);
    end

n = size(sequence, 2);

E = 0;
for i = 1:n-1
    E = E + metric(cities(sequence(i),:), cities(sequence(i+1),:));
end

E = E + metric(cities(sequence(end),:), cities(sequence(1),:));

end
