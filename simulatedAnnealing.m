function [result] = simulatedAnnealing(cities, initialTemperature, endTemperature)
n = size(cities,1); % получаем размер вектора городов

state = randperm(n); % задаём начальное состояние, как случайный маршрут
% Функция randperm(n) - генерирует случайныую последовательность из целых чисел от 1 до n
delete(findall(0,'Type','figure'))
scatter(cities(:, 1), cities(:, 2));
hold on;
for i = 1:n-1
    line([cities(state(i), 1), cities(state(i+1), 1)], [cities(state(i), 2), cities(state(i+1), 2)]);
end
line([cities(state(1), 1), cities(state(end), 1)], [cities(state(1), 2), cities(state(end), 2)]);
hold off;

currentEnergy = calculateEnergy(state, cities); % вычисляем энергию для первого состояния
currentEnergy
T = initialTemperature;

for i = 1:10000  %на всякий случай ограничеваем количество итераций
    % может быть полезно при тестировании сложных функций изменения температуры T 
    stateCandidate = generateStateCandidate(state.').'; % получаем состояние-кандидат
    candidateEnergy = calculateEnergy(stateCandidate, cities); % вычисляем его энергию
    
    if(candidateEnergy < currentEnergy) % если кандидат обладает меньшей энергией
        currentEnergy = candidateEnergy; % то оно становится текущим состоянием
        state = stateCandidate;
    else
        p = getTransitionProbability(candidateEnergy-currentEnergy, T); % иначе, считаем вероятность
        if (isTransition(p)) % и смотрим, осуществится ли переход
            currentEnergy = candidateEnergy;
            state = stateCandidate;
        end
    end;
    
    T = decreaseTemperature(initialTemperature, i); % уменьшаем температуру
    if(T <= endTemperature) % условие выхода
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
%initialTemperature - начальная температура
%i - номер итерации
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
%seq - предыдущее состояние (маршрут), из которого
%мы хотим получить состояние-кандидат

n = size(seq,1); % определяем размер последовательности
l = randi(n,1); % генерируем целое случайное число
r = randi(n,1);% генерируем целое случайное число
res = seq;
if(l > r)
    res(r:l) = flipud(seq(r:l)); % обращаем подпоследовательность
else
    res(l:r) = flipud(seq(l:r)); % обращаем подпоследовательность
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
