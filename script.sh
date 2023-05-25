#Gabija Genčiūtė
#Parsisiunčiu failus iš ENA
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR204/ERR204044/ERR204044_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR204/ERR204044/ERR204044_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/030/SRR15131330/SRR15131330_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/030/SRR15131330/SRR15131330_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR182/064/SRR18214264/SRR18214264_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR182/064/SRR18214264/SRR18214264_2.fastq.gz
#tada atlieku fastqc
fastqc inputs/* -o outputs/raw_fastqc
#ERR204044_2 readas pagal fastqc rodo šauktuką(šiek tiek blogesnė, nei "tobula" kokybė), kiti readai iš esmės yra labai geros kokybės, bet pradžios ir pabaigos
#galima sakyti, kad yra blogesnės kokybės, todėl spėju, kad po kirpimo, bus geresnė kokybė :) labai blogos kokybės readų - nėra.
#tada atlieku read'ų kirpimą
thrds=6
for i in inputs/*_1.fastq.gz;
do 
    R1=${i}
    R2="inputs/"$(basename $i _1.fastq.gz)"_2.fastq.gz"
    TO="outputs/trimmed/"
    trim_galore -j ${thrds} -o ${TO} --paired ${R1} ${R2};
done
fastqc outputs/trimmed/*val* outputs/val_fastqc
#spėjimas pasitvirtino - po kirpimo SRR151* readai tapo galima sakyti tobuli, problemiškas yra ERR204044_2 readas (vis dar pradžioje yra vidutinė kokybė)
#bet viskas iš esmės yra gerai, o SRR182* readams fastqc rodo žalią varnelę - nors pabaigoje kokybė truputį peržengia mažiau 28 ribą.
multiqc outputs/ .
#multiqc report'ą įdedu į git'ą.
#############################################################################
#vėliau naudosiu gal
VR1="outputs/trimmed/"$(basename ${i} _1.fastq.gz)"_1_val_1.fq.gz"
VR2="outputs/trimmed/"$(basename ${i} _1.fastq.gz)"_2_val_2.fq.gz"


