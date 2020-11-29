from scriptutils import *
import pprint
#Read go annotations
gotext = open("GOannotations.txt", "r")
goannotations = []
next(gotext)
for line in gotext:
    content_array = line.split(maxsplit=2)
    annotation_element = {
        "gene_stab_id":content_array[0],
        "go_term":content_array[1],
        "term_name":get_list_elem(content_array, 2)
    }
    goannotations.append(annotation_element)
gotext.close()

#Read all genes
allgenes = create_list("allGenes(2).txt")

#Read ids2
ids2genes = create_list("IDS2.txt")

#For generating annotations
"""goids = []
for gene in goannotations:
    goids.append(gene['gene_stab_id'])

f = open('annotadedgenes.csv', 'w')
f.write('id,go,sel\n')

for gene in allgenes:
    go = gene in goids
    sel = gene in ids2genes
    f.write(f'{gene},{go},{sel}\n')
    print(f'{gene},{go},{sel}\n')
"""

contingency_table = [[0,0],[0,0]]
file_annotated = open('annotadedgenes.csv', 'r')
next(file_annotated)
for line in file_annotated:
    go = line.split(',')
    if go[1] == 'True' and go[2].rstrip() == 'True':
        contingency_table[0][0] += 1
    elif go[1] == 'True' and go[2].rstrip() == 'False':
        contingency_table[0][1] += 1
    elif go[1] == 'False' and go[2].rstrip() == 'True':
        contingency_table[1][0] += 1
    elif go[1] == 'False' and go[2].rstrip() == 'False':
        contingency_table[1][1] += 1
file_annotated.close()

expected = calculate_chi_expected_values(contingency_table)
print(calculate_chi_square(contingency_table,expected))
