Olá!
Para essa automação funcionar, uma série de critérios precisa atingido:


Primeiro: O poconv, psrmin,psh0 e psh1 precisam estar compilados e funcionando dentro da pasta arquivos

Segundo: o MSCD-ATA também precisa estar compilado e funcionando dentro da pasta arquivos

Terceiro: o arquivo de input precisa sempre ter o nome "input_cluster.txt" ele é a referência para criar o cluster.i

Quarto: o arquivo experimental precisa estar dentro da pasta "arquivos" e estar no modelo do MSCD

Quinto: para o new.py que plota o padrão teórico e cálcula o R-factor basta deixar o nome do arquivo de saída de saida.out.
Ou entre dentro do new.py e mude lá nas últimas linhas o nome do arquivo q ele vai ler.

Sexto: Use o lattice como "1.000" no cabeçalho do "input_cluster" dessa forma os Ps e Rm serão cálculados corretamente

Sétimo: Por limitação do programa q cálcula a diferença de fase, o input_cluster.txt só pode ter uma base de vetores, mas nada impede vc de usar essa automação apenas para criar o os ps e rm e depois pegar o arquivo de input para usar nas suas simulações

Para plotar a estrutura, execute python3 plot_structure.py em seguida dê VESTA no terminal e abra o arquivo cluster.xyz

