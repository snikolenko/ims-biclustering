alpha = 0.5; %акцент пространственных ребер в матрице смежности
beta = 1 - alpha; %акцент двудольных ребер в матрице смежности
max_force = 120; %вес ребра на соседних точках
otsechka = 1.5; %расстояние, в пределах которого вводятся пространственные ребра 
k = round(8*spectors/9);
r = 0.20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% на всякий случай все очищаем
clear specs_diag;
clear pix_diag;
clear W;
clear bigW;
clear bigD;
clear bigL;
clear V;
clear maximums;






%%%%%%%%%%%%%%%%%%%%%%%%%%% вычисляем все вспомогательные величины %%%%%%%%

maximums(1:pixels) = 0;
for i= 1:pixels
    tmp=0;
    for j=1:spectors
        if spectra(j, i) > tmp
           tmp = spectra(j, i);
        end
    end
    maximums(i) = tmp;
end %в maximums храним максимумы яркости для каждой точки

specs_diag(1:spectors)=0;
pix_diag(1:pixels)=0;
for i=1:spectors
    for j=1:pixels
        specs_diag(i)=specs_diag(i) + spectra(i,j); %тут храним суммы всех двудольных ребер для m/z-значений
        pix_diag(j)=pix_diag(j) + spectra(i,j); %тут храним суммы всех двудольных ребер для точек
    end
end

W(1:pixels, 1:pixels) = 0; % квадратная матрица с пространственными весами
on_diag(1:pixels) = 0; % 
little(1:pixels) = 0; % 

Z(1:pixels, 1:spectors) = 0;
I(1:pixels, 1:spectors) = 0;
U(1:pixels, 1:spectors) = 0;

for i=1:pixels
    [Z(i, :), I(i, :)] = sort(spectra(:, i));
    for j = k+1:spectors
        U(i, I(i, j)) = 1;
    end
end %в U будем хранить в бинарной форме те m/z-спектры, что лежат в списке сильных для каждой точки

clear Z;
clear I;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%%%%%%%%%%%%%%%%%%% заполняем матрицу смежности %%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:pixels
    x_i = x_printed(i);
    y_i = y_printed(i);
    for j=1:i-1
        x_j = x_printed(j);
        y_j = y_printed(j);
%         
        dist = sqrt((x_i - x_j)^2 + (y_i - y_j)^2); % расстояние
        x_comp = sign(x_i - x_j);
        y_comp = sign(y_i - y_j);
        if dist <= otsechka % если оно меньше отсечки
            dt = dot(U(i, :), U(j, :)) / (spectors - k); % мера похожести пикселей - доля общих ярких m/z значений
            W(i,j) = alpha * max_force * dt / dist; % вес пространственного ребра пропорционален похожести точек и обратно пропорционален расстоянию
            
            W(j,i) = W(i,j);
            on_diag(i) = on_diag(i) + W(i,j);
            on_diag(j) = on_diag(j) + W(j,i);
            if dt <= r
                little(i) = little(i) + W(i,j);
                little(j) = little(j) + W(j,i);
            end
        end
        
    end
end



bigD(1:(pixels+spectors), 1:(pixels+spectors))=0; % диагональная матрица
bigW(1:(pixels+spectors), 1:(pixels+spectors))=0; % полная матрица смежности
for i=1:pixels
    bigD(i,i)=on_diag(i) + beta*pix_diag(i);
    bigW(i,i) =little(i) +  beta*pix_diag(i);
end
for i=1:spectors
    bigD(pixels + i,pixels + i)=beta*specs_diag(i);
    bigW(pixels + i,pixels + i)=beta*specs_diag(i);
end

bigM=[W beta*transpose(spectra); beta*spectra zeros(spectors)]; % полная матрица смежности всего графа
clear W;

bigL = bigD - bigM; % вычислили матрицу Лапласа графа
clear bigM;
clear bigD;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[V, D] = eigs(bigL, bigW, 10, 'sa'); % в V - первые 10 собственных векторов, в D - их собственные значения

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%%%%% визуализация %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear tmp;
clear ID;
m = 10;
vecs = 7;
tmp=V(1:spectors + pixels, 2:vecs+1);
ID = kmeans(tmp,m,...
                'Options', statset('MaxIter', 250));

size = 20;
pic(1:pixels, 3)=0;
for i=1:pixels
    pic(i,1)=x_printed(i);
    pic(i,2)=y_printed(i);
    pic(i,3)=ID(i);
end

figure
G=scatter(pic(:,2), pic(:,1), size, pic(:,3), 'square', 'fill');

