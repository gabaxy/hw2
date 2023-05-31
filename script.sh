#!/bin/sh
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
#genomų surinkimas naudojant spades ir abyss
gunzip outputs/trimmed/*.gz
mkdir ERR204044_abyss SRR15131330_abyss SRR18214264_abyss
for i in inputs/*_1.fastq.gz;
do 
    thrds=6
    VR1="outputs/trimmed/"$(basename ${i} _1.fastq.gz)"_1_val_1.fq"
    VR2="outputs/trimmed/"$(basename ${i} _1.fastq.gz)"_2_val_2.fq"
    TO="outputs/assembly/"
    spades -t ${thrds} -1 ${VR1} -2 ${VR2} -o ${TO}/"$(basename ${i} _1.fastq.gz)_spades"
    #abyss-pe k=64 name=$(basename ${i} _1.fastq.gz) -C ${TO}$(basename ${i} _1.fastq.gz)_abyss in=${VR1} ${VR2}
    #nelabai man pavyko su šita eilute, todėl atskirai padariau kiekvieną
done

abyss-pe k=64 name=ERR204044 -C ~/HW2/outputs/assembly/ERR204044_abyss in='~/HW2/outputs/trimmed/ERR204044_1_val_1.fq ~/HW2/outputs/trimmed/ERR204044_2_val_2.fq'
abyss-pe k=64 name=SRR15131330 -C ~/HW2/outputs/assembly/SRR15131330_abyss in='~/HW2/outputs/trimmed/SRR15131330_1_val_1.fq ~/HW2/outputs/trimmed/SRR15131330_2_val_2.fq'
abyss-pe k=64 name=SRR18214264 -C ~/HW2/outputs/assembly/SRR18214264_abyss in='~/HW2/outputs/trimmed/SRR18214264_1_val_1.fq ~/HW2/outputs/trimmed/SRR18214264_2_val_2.fq'

#tada parsisunčiu ref. seką CP015498 ir paleidžiu per usegalaxy.eu QUAST.
#QUAST parodė, jog patikimiausios sekos(pagal pateiktą heatmap'ą - mėlyniausios) yra sekos, kurios buvo surinktos SPADES programa. Šios sekos buvo didesniu procentu padengusios ref. genomą(77-82%)
#Lyginant su abyss programa, kur rezultatai buvo 73-76%. Taip pat su spades programa rinktų genomų NG yra didesnis(geras ženklas), ir LG yra mažesnis(geras ženklas)
#SRR15131330 read'as turėjo panašius rezultatus ir renkant abyss programa, tačiau padengtumas 5% mažesnis reiškė, jog visi kiti rezultatai yra taip pat mažesni -
#ilgis, kuris buvo sulygintas su ref. genomu, ir kiti rezultatai(turiu omenyje, kad tai yra logiška)
#blogiausius rezultatus pasiekė SRR18214264, rinktas abyss programa - pats mažiausias procentas - beveik 74%, sulygintas su ref. genomu ilgis - 1509876, darant su spades - 1534348
#reziume - Spades programa sulygino daug geriau, nei abyss, nors ir kai kurie rezultatai abyss programa darant buvo šiek tiek geresni. Todėl map'insiu ir toliau naudosiu tik spades failus

for i in inputs/*_1.fastq.gz;
do
    thrds=6
    VR1="outputs/trimmed/"$(basename ${i} _1.fastq.gz)"_1_val_1.fq"
    VR2="outputs/trimmed/"$(basename ${i} _1.fastq.gz)"_2_val_2.fq"
    CTG="outputs/assembly/"$(basename ${i} _1.fastq.gz)"_spades/contigs.fasta"
    bwa index ${CTG}
    bwa mem -t ${thrds} ${CTG} ${VR1} ${VR2} | samtools view -b -@ ${thrds} - > outputs/map/$(basename ${i} _1.fastq.gz).bam
    samtools sort -@ ${thrds} outputs/map/$(basename ${i} _1.fastq.gz).bam -o outputs/map/$(basename ${i} _1.fastq.gz).sorted.bam
    samtools index outputs/map/$(basename ${i} _1.fastq.gz).sorted.bam
    samtools flagstat outputs/map/$(basename ${i} _1.fastq.gz).sorted.bam > outputs/map/$(basename ${i} _1_val_1.fq.gz)_map_stats.txt
done
#fractions ir coverage(naudojant awk įrankį)
#ERR204044 read'as: 99,62% dalies, 295.777
#SRR15131330 read'as: 88.82% dalies, 1836.51
#SRR18214264 read'as: 99.73% dalies, 276,785
#ERR204* ir SRR182* - labai geri rezultatai, SRR151 - vidutiniškai, o coverage - labai per didelis??


#Gepard įrankis man sukūrė 3 nuotraukas, kurias pridedu užduotyje. jose atvaizduojamas lyginimas mano genomų vienas su kitu. geriausias sulyginimas yra tarp ERR204* ir SRR182
#palyginimo nuotraukoje, linija prasideda (0,0) koordinatėje, tai reiškia, kad tai yra tiesioginis atitikimas. Kiti du palyginimai buvo mažiau taiklūs atitikimo atžvilgiu, nes linija neprasidėjo (0,0) taške
#Padarytas RINGPLOT(iš tiesų - nežinau, ar teisingas, nes lyginau su pateiktu ref. genomu rast duombazėj) rodo identiškus rezultatus, kaip quast - SRR182* atitikimas ir padengtumas yra didesnis
#nei likusių dviejų, o likę du - panašesni tarpusavyje - prastesni rezultatai, tačiau praktiškai identiški vienas kitam.
#Šalia gepard ir ringplot taip pat pridedu nuotrauką, kurioje matosi SRR182* ir ERR204* panašumas jau QUAST programoje plot'uose, kai linijos eina iš esmės ta pačia trajektorija.
#
#nusikopijuoju naudojamus spades genomus į kitą direktoriją, kad galėčiau sukurti duombazes
mkdir outputs/DB
for i in inputs/*_1.fastq.gz;
do
    CTG="outputs/assembly/"$(basename ${i} _1.fastq.gz)"_spades/contigs.fasta"
    cp ${CTG} outputs/DB/"$(basename ${i} _1.fastq.gz)"-contigs.fasta
    makeblastdb -in outputs/DB/$(basename ${i} _1.fastq.gz)-contigs.fasta -dbtype nucl -parse_seqids;
    makeblastdb -in outputs/DB/$(basename ${i} _1.fastq.gz)-contigs.fasta -dbtype prot -parse_seqids;
    blastn -query refs/GCF_022832545.1_ASM2283254v1_genomic.fna -db outputs/DB/$(basename ${i} _1.fastq.gz)-contigs.fasta > outputs/DB/$(basename ${i} _1.fastq.gz)-blastgene_info.txt
    tblastn -query refs/GCF_022832545.1_ASM2283254v1_protein.faa -db outputs/DB/$(basename ${i} _1.fastq.gz)-contigs.fasta > outputs/DB/$(basename ${i} _1.fastq.gz)-blastprot_info.txt
done
#tuomet skaičiuoju nuspėjamų genų kiekį ir overlap. RAST kadangi yra naršyklėje, tai iš ten tiesiog paimu nurodytus skaičius, GeneMarks failai mano pašte,
#todėl apskaičiuoju naudodama ctrl+f('gene_id') ir išveda kiek iš viso tokių eilučių yra - tai lygu kiek yra genų, o blast rezultatus apskaičiuoju linux aplikoje, naudodama awk '{ print $1 }' blastgene_info.txt | sort -u | wc -l
#genemarks: ERR204* 2310 genai, SRR151* 2552 genai, SRR182* 2364 genas rastas
#RAST: ERR204* 2494 genai, SRR151* 2700 genai, SRR182* 2451 genai rasti
#BLAST: ERR204* 1873 genai, SRR151* 2380 genai, SRR182* 1500 genai rasti

#BUSCO:
#SRR264 99 proc. complete ir single-copy sekų, 0.5 proc. fragmentuotų, 0.5 proc. trūkstamų.
#SRR1330 98.5 proc. complete ir single-copy sekų, 0.5 proc. fragmentuotų, 1 proc. trūkstamų.
#ERR 99 proc. complete ir single-copy sekų, 0.5 proc. fragmentuotų, 0.5 proc. trūkstamų.
#apskaičiuoti procentai iš summary nuotraukų, kurias pridedu taip pat į užduotį pavadinimais:busco_*

#taigi išvada: SRR182 ir ERR204 yra labai panašiai surinkti genomai tarpusavyje pagal gepard, quest ir ringplot.

#CITATIONS:

#The RAST Server: Rapid Annotations using Subsystems Technology.
#Aziz RK, Bartels D, Best AA, DeJongh M, Disz T, Edwards RA, Formsma K, Gerdes S, Glass EM, Kubal M, Meyer F, Olsen GJ, Olson R, Osterman AL, Overbeek RA, McNeil LK, Paarmann D, Paczian T, Parrello B, Pusch GD, Reich C, Stevens R, Vassieva O, Vonstein V, Wilke A, Zagnitko O.
#BMC Genomics, 2008, [ PubMed entry ]
#The SEED and the Rapid Annotation of microbial genomes using Subsystems Technology (RAST).
#Overbeek R, Olson R, Pusch GD, Olsen GJ, Davis JJ, Disz T, Edwards RA, Gerdes S, Parrello B, Shukla M, Vonstein V, Wattam AR, Xia F, Stevens R.
#Nucleic Acids Res. 2014 [ PubMed entry ]
#RASTtk: A modular and extensible implementation of the RAST algorithm for building custom annotation pipelines and annotating batches of genomes.
#Brettin T, Davis JJ, Disz T, Edwards RA, Gerdes S, Olsen GJ, Olson R, Overbeek R, Parrello B, Pusch GD, Shukla M, Thomason JA, Stevens R, Vonstein V, Wattam AR, Xia F.
#Sci Rep., 2015, [ PubMed entry ]

#usegalaxy.eu:
#Quast
#- Mikheenko, A., Prjibelski, A., Saveliev, V., Antipov, D., & Gurevich, A. (2018). Versatile genome assembly evaluation with QUAST-LG. Bioinformatics, 34(13), i142–i150. https://doi.org/10.1093/bioinformatics/bty266  
#- Mikheenko, A., Valin, G., Prjibelski, A., Saveliev, V., & Gurevich, A. (2016). Icarus: visualizer for de novo assembly evaluation. Bioinformatics, 32(21), 3321–3323. https://doi.org/10.1093/bioinformatics/btw379  
#- Gurevich, A., Saveliev, V., Vyahhi, N., & Tesler, G. (2013). QUAST: quality assessment tool for genome assemblies. Bioinformatics, 29(8), 1072–1075. https://doi.org/10.1093/bioinformatics/btt086  
#RagTag
#- Alonge, M., Soyk, S., Ramakrishnan, S., Wang, X., Goodwin, S., Sedlazeck, F. J., Lippman, Z. B., & Schatz, M. C. (2019). RaGOO: fast and accurate reference-guided scaffolding of draft genomes. 20(1). https://doi.org/10.1186/s13059-019-1829-6  
#BUSCO
#- Simão, F. A., Waterhouse, R. M., Ioannidis, P., Kriventseva, E. V., & Zdobnov, E. M. (2015). BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. Bioinformatics, 31(19), 3210–3212. https://doi.org/10.1093/bioinformatics/btv351 

#GeneMarkS-2
#Lomsadze A, Gemayel K, Tang S, Borodovsky M
#Modeling leaderless transcription and atypical genes results in more accurate gene prediction in prokaryotes.
#Genome Res, 2018, 29(7), pp 1079-1089


