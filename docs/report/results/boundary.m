
projroot = '../../../';
files = [projroot  'build/analysis2/'];
target = 'gen2/';


graphics_toolkit('gnuplot')

mbrot = dlmread([files 'mbrot0.dat'], ' ', 1, 0);

for n = 1:11
    mbrotp = dlmread([files 'mbrot' num2str(n) '.dat'], ' ', 1, 0);

    mbrot = mbrot+mbrotp;
end

mbrot = mbrot ./ 12;

dlmwrite([target 'mbrot-total.dat'], mbrot, ' ');
