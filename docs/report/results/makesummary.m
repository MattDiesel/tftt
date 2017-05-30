
projroot = '../../../';
files = [projroot  'build/analysis/'];
target = 'gen/';


graphics_toolkit('gnuplot')


mbrotFiles = 0:11;
worldSizes = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,20,22,24,26,28,30,32];

summry = zeros(length(worldSizes), length(mbrotFiles)+1);

summry(:,1) = worldSizes';


for n = mbrotFiles
    mbrot = dlmread([files 'mbrot' num2str(n) '.dat'], ' ', 1, 0);

    ccells = mbrot(11,1) ** 0.4562;

    for w = 1:length(worldSizes)
        summry(w,n+2) = mbrot(11,w+1) / ccells;
    end
end

dlmwrite([target 'summ2.dat'], summry, ' ');
