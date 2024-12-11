function I = random_generator (PDF)

%% Normalized check

PDF = PDF ./ sum(PDF,'all');
PDF_dimension = size(PDF);

%%

PDF_linear = reshape(PDF,numel(PDF),1);
PDF_linear = cumsum(PDF_linear./sum(PDF_linear));

random_number = rand;
random_index = find (random_number < PDF_linear,1,'first');

[I{1:numel(PDF_dimension)}] = ind2sub(PDF_dimension,random_index);
I = cell2mat(I)-1;