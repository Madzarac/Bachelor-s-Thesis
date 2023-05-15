import re


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

def get_clipping_length(cigar_string):

    clipping = 0
    pattern = re.compile(r'(\d+)([MIDNSHPX=])')
    cigar_tuples = pattern.findall(cigar_string)

    if cigar_tuples[0][1] in ('S', 'H'):
        clipping += int(cigar_tuples[0][0])

    if cigar_tuples[-1][1] in ('S', 'H'):
        clipping += int(cigar_tuples[-1][0])

    return (clipping)

def analyse(lines, target_rank, taxonomy_tree, ranks):
    results = {}
    max_values = {}
    transformed_rows = {}
    counter4 = 0
    for line in lines:
        parts = re.split(r'\t+', line.strip())
        if len(parts) < 11:
            continue
        if parts[1].strip() == "4": # int(parts[4].strip()) <= 45: #FLAG 4 means read is unmapped, we're also not going to consider reads with MAPQ value less than 30
            counter4 += 1
            continue
        nmPart = re.split(r':+', parts[11].strip())
        nm = nmPart[-1]
        if get_clipping_length(parts[5]) > 100 and int(parts[4].strip()) < 45 and int(nm) > len(parts[9]) * 0.001:
            continue
        read_id = parts[0].strip()
        tax_id_extended = parts[2].strip()

        parts2 = re.split(r'\|+', tax_id_extended.strip())
        if len(parts2) == 1:
            continue
        tax_id = parts2[2].strip()
        resulting_tax_id = find_resulting_tax_id(tax_id, target_rank, taxonomy_tree, ranks)
        if resulting_tax_id == 0:
            resulting_tax_id = tax_id

        value_cig = 0.0
        if len(parts) >= 14 :           #AS filed: alingment score 
            aS = parts[13].strip()      #NM field: Edit distance to the reference (number of mismatches), does not count clippings
            parts3 = re.split(r':+', aS.strip())
            value_cig = int(parts3[-1].strip()) - int(nm) - get_clipping_length(parts[5]) + int(parts[4].strip())


        if read_id in results:
            if value_cig > max_values[read_id]:
                results[read_id] = []
                results[read_id].append(resulting_tax_id)
                max_values[read_id] = value_cig
            elif value_cig == max_values[read_id]:
                results[read_id].append(resulting_tax_id)
        else:
            results[read_id] = []
            results[read_id].append(resulting_tax_id)
            max_values[read_id] = value_cig

    for read_id in results:
        tax_ids = results[read_id]
        values = {}
        for tax_id in tax_ids:
            if tax_id not in values:
                values[tax_id] = 0
            values[tax_id] += 1
        max_tax_id = ""
        max_value = 0
        for tax_id in values:
            value = values[tax_id]
            if value > max_value:
                max_value = value
                max_tax_id = tax_id
        transformed_rows[read_id.strip()] = (max_tax_id, ranks[max_tax_id])
    print(counter4)
    return transformed_rows


def main_func(number_of_nodes_files, path):
    ranks_lists = {}
    taxonomy_tree_lists = {}
    taxonomy_names_lists = {}

    for nodes_index in range(0, int(number_of_nodes_files)):
        nodes_file = open(path + "/" + "nodes" + str(nodes_index) + ".dmp", "r")
        nodes_lines = nodes_file.readlines()
        ranks = {}
        taxonomy_tree = {}
        for line in nodes_lines:
            parts = re.split(r'\t\|\t+', line.strip())
            taxonomy_tree[parts[0].strip()] = parts[1].strip()
            ranks[parts[0].strip()] = parts[2].strip()
        ranks_lists[nodes_index] = ranks
        taxonomy_tree_lists[nodes_index] = taxonomy_tree

    for names_index in range(0, int(number_of_nodes_files)):
        names_file = open(path + "/" + "names" + str(names_index) + ".dmp", "r")
        names_lines = names_file.readlines()
        taxonomy_names = {}
        for line in names_lines:
            parts = re.split(r'\|+', line.strip())
            if parts[3].strip() == "scientific name":
                taxonomy_names[parts[0].strip()] = parts[1].strip()
        taxonomy_names_lists[names_index] = taxonomy_names


    target_ranks = ["species", "genus"]
    results = "new_results"

    databases = ["velika_baza_bakterija", "mala_baza1", "mala_baza2"]
    datasets = [
        ("zymo_ont.badread.fastq.gz", 1),
        ("zymo_ont_veliki.badread.fastq.gz", 2)
    ]


    for (dataset, number_of_dataset) in datasets:
        print("#Dataset: " + str(dataset))

        for database in databases:
            if number_of_dataset == 1 and database == "mala_baza2":
                continue
            if number_of_dataset == 2 and database == "mala_baza1":
                continue
            print("#Database: " + str(database))
            for target_rank in target_ranks:
                print("Target rank: " + str(target_rank))

                results_filename = "aln" + str(number_of_dataset) + ".sam"
                results_file = open(results_filename, "r")
                results_lines = results_file.readlines()

                parsed_rows = []
                taxonomy_tree = taxonomy_tree_lists[0]
                ranks = ranks_lists[0]

                parsed_rows = analyse(results_lines, target_rank, taxonomy_tree, ranks)

                transformed_rows = parsed_rows
                filename = results + "/" + str(database) + "_" + str(number_of_dataset) + "_" + str(target_rank) + ".f2"
                outfile = open(filename, "w")

                for read_id in transformed_rows:
                            (tax_id, rank) = transformed_rows[read_id]
                            outfile.write(read_id.strip() + "\t" + tax_id.strip() + "\t" + rank.strip() + "\n")

    results_file.close()
    outfile.close()



main_func(3, "names_nodes")
