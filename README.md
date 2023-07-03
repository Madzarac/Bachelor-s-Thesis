# Bachelor-s-Thesis
Klasifikacija metagenomskog uzorka pomoću alata za mapiranje
Metagenomika je relativno nova grana bioinformatike koja se bavi proučavanjem 
genetskih uzoraka prikupljenih direktno iz okoliša te nam nudi mogućnost primjene 
stečenog znanja u mnogim područjima, među kojima su medicina, agrikultura i ekologija. 
Jedan od osnovnih koraka u analizi metagenomskog uzorka je određivanje organizama koji 
se u njemu nalaze te njihove zastupljenosti u uzorku. Iz toga razloga važno je koristiti 
dobre alate za klasifikaciju. Ovaj rad bavi se osmišljavanjem algoritma koji na temelju
mapiranja pokušava odrediti sastav uzorka. Pomoću alata Kraken2 izgradile su se dvije 
baze bakterija te se na njih mapiraju očitanja korištenjem alata minimap2. Nad dobivenim 
mapiranjima provela se klasifikacija pomoću osmišljenog algoritma te su dobiveni rezultati 
uspoređeni s rezultatima klasifikacije dobivenim uporabom alata Kraken2.

Metagenomic sample classification using mapping tools
Metagenomics is a relatively new branch of bioinformatics that focuses on studying 
genetic samples collected directly from the environment. It has potential to be used to 
improve various fields, such as medicine, agriculture and ecology. One of the fundamental
steps in analyzing a metagenomic samlple is determining the organisms present in the 
sample. Therefore, it is important to use reliable classification tools. This paper involves 
the development o fan algorithm that utilizes mapping to determine the composition of a 
sample. First, two custom bacterial databases we constructed using the Kraken 2 tool, and 
the reads were mapped onto these databases using the tool Minimap 2. Classification was 
performed based on the generated mappings using the developed algorithm, and the results 
were compared with the classification results obtained using the Kraken2 tool.
