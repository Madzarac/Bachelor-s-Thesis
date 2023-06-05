import re
import sys

def find_resulting_tax_id(current_tax, target_rank, taxonomy_tree, ranks):
    not_found_resulting_tax_id = True
    resulting_tax_id = 0
    tax = current_tax
    while not_found_resulting_tax_id:
        if tax in taxonomy_tree:
            parent_rank = ranks[tax]
            if parent_rank == target_rank:
                not_found_resulting_tax_id = False
                resulting_tax_id = tax
            else:
                prev_tax = tax
                tax = taxonomy_tree[tax]
                if tax == prev_tax:
                    not_found_resulting_tax_id = False
        else:
            not_found_resulting_tax_id = False
    return resulting_tax_id


def analyse(lines, target_rank, taxonomy_tree, ranks, isSam):
    results = {}
    max_values = {}
    transformed_rows = {}
    for line in lines:
        resulting_tax_id = 0
        value = 0.0
        parts = re.split(r'\t+', line.strip())
        read_id = parts[0].strip()
        if line.startswith("@"):
            continue
        if isSam:
            if parts[1].strip() == "4": 
                continue
            
            tax_id_extended = parts[2].strip()

            parts2 = re.split(r'\|+', tax_id_extended.strip())
            if len(parts2) == 1:
                continue
            tax_id = parts2[2].strip()
            resulting_tax_id = find_resulting_tax_id(tax_id, target_rank, taxonomy_tree, ranks)
            if resulting_tax_id == 0:
                resulting_tax_id = tax_id
            
            if len(parts) >= 11 :          # making sure that the opptional fields are included  
                nmPart = re.split(r':+', parts[11].strip())  # NM field: Edit distance to the reference (number of mismatches), does not count clippings
                nm = nmPart[-1]
                value = int(nm)

        else: #for .paf files
            nmPart = re.split(r':+', parts[12].strip())
            nm = nmPart[-1]
            value = int(nm)

            tax_id_extended = parts[5].strip()
            parts2 = re.split(r'\|+', tax_id_extended.strip())
            if len(parts2) == 1:
                continue
            tax_id = parts2[2].strip()
            resulting_tax_id = find_resulting_tax_id(tax_id, target_rank, taxonomy_tree, ranks)
            if resulting_tax_id == 0:
                resulting_tax_id = tax_id

        # if one read has multiple alignments, select one with the greatest assigned value
        if read_id in results:
            if value < max_values[read_id]:
                results[read_id] = resulting_tax_id
                max_values[read_id] = value
        else:
            results[read_id] = resulting_tax_id
            max_values[read_id] = value

    for read_id in results:
        tax_ids = results[read_id]
        transformed_rows[read_id.strip()] = (tax_ids, ranks[tax_ids])

    return transformed_rows


def main_func(database, resultsFile, file):

    nodes_file = open(database + "/taxonomy/nodes.dmp", "r")
    nodes_lines = nodes_file.readlines()
    ranks = {}
    taxonomy_tree = {}
    for line in nodes_lines:
        parts = re.split(r'\t\|\t+', line.strip())
        taxonomy_tree[parts[0].strip()] = parts[1].strip()
        ranks[parts[0].strip()] = parts[2].strip()

    names_file = open(database + "/taxonomy/names.dmp", "r")
    names_lines = names_file.readlines()
    taxonomy_names = {}
    for line in names_lines:
        parts = re.split(r'\|+', line.strip())
        if parts[3].strip() == "scientific name":
            taxonomy_names[parts[0].strip()] = parts[1].strip()


    target_ranks = ["species", "genus"]

    for target_rank in target_ranks:
        print("Target rank: " + str(target_rank))

        results_file = open(file, "r")
        results_lines = results_file.readlines()

        parsed_rows = []
                
        isSam = False
        if file.endswith(".sam"):
            isSam = True
        parsed_rows = analyse(results_lines, target_rank, taxonomy_tree, ranks, isSam)

        transformed_rows = parsed_rows
        filename = resultsFile + "_" + str(target_rank) + ".f2"
        outfile = open(filename, "w")

        for read_id in transformed_rows:
            (tax_id, rank) = transformed_rows[read_id]
            outfile.write(read_id.strip() + "\t" + tax_id.strip() + "\t" + rank.strip() + "\n")

    results_file.close()
    outfile.close()


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("make sure to provide all neded arguments: database, where do you want the results, file in sam or paf format.")
    if not (sys.argv[3].endswith(".sam") or sys.argv[1].endswith(".paf")):
        print("File has to be in SAM or PAF format.")
    else:
        main_func(sys.argv[1], sys.argv[2], sys.argv[3])
