import re
import pysam

folder = "krakenPosNegResults"
realTaxidsFile = []
realTaxidsFile.append("match_id_taxid/mala_baza1_taxids")
realTaxidsFile.append("match_id_taxid/mala_baza2_taxids")
krakenSequences = []
krakenSequences.append("sequences1.kraken")
krakenSequences.append("sequences2.kraken")



for i in range(0,2):
    dictKraken = {}
    dictKrakenUnclass = {}
    dictRealTaxidValues = {}

    falsePositive = {}
    falseNegative = {}
    truePositive = {}
    trueNegative = {}

    real = open(realTaxidsFile[i], 'r')
    realTaxids = real.readlines()

    species = open(krakenSequences[i], 'r')
    speaciesTaxid = species.readlines()

    falsePos = open(folder + "/" + "falsePositive" + str(i + 1) + ".txt", "w")
    falseNeg = open(folder + "/" + "falseNegative" + str(i + 1) + ".txt", "w")
    truePos = open(folder + "/" + "truePositive" + str(i + 1) + ".txt", "w")
    trueNeg = open(folder + "/" + "trueNegative" + str(i + 1) + ".txt", "w")

    for line in realTaxids:
        parts = re.split(r'\t+', line.strip())
        dictRealTaxidValues[parts[0]] = parts[1]

    for line in speaciesTaxid:
        parts = re.split(r'\t+', line.strip())
        if parts[0] == "U":
            taxid = parts[2]
            dictKrakenUnclass[parts[1]] = taxid
        else:
            taxid = parts[2]
            dictKraken[parts[1]] = taxid


    correct = 0
    for id in dictKraken.keys():
        if id in dictRealTaxidValues:
            if dictRealTaxidValues[id] == dictKraken[id]:
                correct += 1
                truePositive[id] = "\t" + dictRealTaxidValues[id] + "\t" + dictKraken[id]
            else:
                falsePositive[id] = "\t" + dictRealTaxidValues[id] + "\t" + dictKraken[id]
        else:
            falsePositive[id] = "\t-\t" + dictKraken[id]

    for id in dictKrakenUnclass.keys():
        if dictRealTaxidValues[id] == "?":
            trueNegative[id] = "\t-" + "\t*"
        else:
            falseNegative[id] = "\t" + dictRealTaxidValues[id] + "\t*"

    falsePos.write("Number of false positive reads: " + str(len(falsePositive)) + ".\n\n")
    falseNeg.write("Number of false negative reads: " + str(len(falseNegative)) + ".\n\n")
    truePos.write("Number of true positive reads: " + str(len(truePositive)) + ".\n\n")
    trueNeg.write("Number of true negative reads: " + str(len(trueNegative)) + ".\n\n")

    falsePos.write("id\tcorrect_species_taxid\tkreken_species_taxid\n")
    falseNeg.write("id\tcorrect_species_taxid\tkreken_species_taxid\n")
    truePos.write("id\tcorrect_species_taxid\tkraken_species_taxid\n")
    trueNeg.write("id\tcorrect_species_taxid\tkreken_species_taxid\n")

    for key, value in falsePositive.items():
        falsePos.write(key + value + "\n")

    for key, value in falseNegative.items():
        falseNeg.write(key + value + "\n")

    for key, value in truePositive.items():
        truePos.write(key + value + "\n")

    for key, value in trueNegative.items():
        trueNeg.write(key + value + "\n")


    falsePos.close()
    falseNeg.close()
    truePos.close()
    trueNeg.close()


