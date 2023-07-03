import re
import sys
from Bio import SeqIO
import pysam


def compare(speciesFile, genusFile, realTaxidsFile, samfile, folder):
    dictSpicies = {}
    dictGenus = {}
    dictRealTaxidValues = {}
    ids = []

    falsePositive = {}
    falseNegative = {}
    truePositive = {}
    trueNegative = {}

    species = open(speciesFile, 'r')
    speaciesTaxid = species.readlines()

    genus = open(genusFile, 'r')
    genusTaxid = genus.readlines()

    real = open(realTaxidsFile, 'r')
    realTaxids = real.readlines()

    digit_index = digit_index = samfile[3]
    falsePos = open(folder + "/" + "falsePositive" + str(digit_index) + ".txt", "w")
    falseNeg = open(folder + "/" + "falseNegative" + str(digit_index) + ".txt", "w")
    truePos = open(folder + "/" + "truePositive" + str(digit_index) + ".txt", "w")
    trueNeg = open(folder + "/" + "trueNegative" + str(digit_index) + ".txt", "w")

    for line in realTaxids:
        sAndGtaxid = []
        parts = re.split(r'\t+', line.strip())
        sAndGtaxid.append(parts[1])
        sAndGtaxid.append(parts[2])
        dictRealTaxidValues[parts[0]] = sAndGtaxid
        ids.append(parts[0])

    for line in speaciesTaxid:
        parts = re.split(r'\t+', line.strip())
        taxid = parts[1]
        dictSpicies[parts[0]] = taxid

    for line in genusTaxid:
        parts = re.split(r'\t+', line.strip())
        taxid = parts[1]
        dictGenus[parts[0]] = taxid

    correctSpecies = 0
    correctGenus = 0
    for id in ids:
        if id not in dictSpicies.keys() and id not in dictGenus.keys():
            if dictRealTaxidValues[id][0] == "?":
                trueNegative[id] = "\t-\t-\t*\t*"
            else: 
                falseNegative[id] = "\t" + dictRealTaxidValues[id][0] + "\t" + dictRealTaxidValues[id][1] + "\t*\t*"
        elif dictRealTaxidValues[id][0] == "?":
            if id not in dictSpicies.keys() and id not in dictGenus.keys():
                trueNegative[id] = "\t-\t-\t*\t*"
            else:
                falsePositive[id] = "\t-\t-" + "\t" + dictSpicies[id] + "\t" + dictGenus[id]
        elif dictRealTaxidValues[id][0] == dictSpicies[id] and dictRealTaxidValues[id][1] == dictGenus[id]:
            truePositive[id] = "\t" + dictRealTaxidValues[id][0] + "\t" + dictRealTaxidValues[id][1] + "\t" + dictSpicies[id] + "\t" + dictGenus[id]
        else: 
            falsePositive[id] = "\t" + dictRealTaxidValues[id][0] + "\t" + dictRealTaxidValues[id][1] + "\t" + dictSpicies[id] + "\t" + dictGenus[id]

    falsePos.write("Number of false positive reads: " + str(len(falsePositive)) + ".\n\n")
    falseNeg.write("Number of false negative reads: " + str(len(falseNegative)) + ".\n\n")
    truePos.write("Number of true positive reads: " + str(len(truePositive)) + ".\n\n")
    trueNeg.write("Number of true negative reads: " + str(len(trueNegative)) + ".\n\n")
    
    falsePos.write("id\tcorrect_species_taxid\tcorrect_genus_taxid\tsam_species_taxid\tsam_genus_taxid\n")
    falseNeg.write("id\tcorrect_species_taxid\tcorrect_genus_taxid\tsam_species_taxid\tsam_genus_taxid\n")
    truePos.write("id\tcorrect_species_taxid\tcorrect_genus_taxid\tscript_species_taxid\tscript_genus_taxid\n")
    trueNeg.write("id\tcorrect_species_taxid\tcorrect_genus_taxid\tscript_species_taxid\tscript_genus_taxid\n")
    
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

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Missing speciesFile, GenusFile, realTaxisFile, samfile and/or folderName")
    else:
        compare(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
