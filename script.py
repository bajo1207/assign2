from scriptutils import *

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
print(contingency_table)


unique_set_of_goterms = set()
for goterm in goannotations:
    unique_set_of_goterms.add(goterm['go_term'])
unique_list_of_goterms = list(unique_set_of_goterms)

print(len(unique_list_of_goterms))

list_of_go_terms_with_chi2 = []
for go_category in unique_list_of_goterms:
    go_ids = []
    contingency_table2 = [[0,0],[0,0]]
    for go in goannotations:
        if go['go_term'] == go_category:
            go_ids.append(go['gene_stab_id'])
    for gene in allgenes:
        go = gene in go_ids
        sel = gene in ids2genes
        if go and sel:
            contingency_table2[0][0] += 1
        elif go and not sel:
            contingency_table2[0][1] += 1
        elif not go and sel:
            contingency_table2[1][0] += 1
        elif not go and not sel:
            contingency_table2[1][1] += 1
    ccot2 = contingency_table2[:]
    expected = calculate_chi_expected_values(ccot2)
    print(expected)
    print(contingency_table2)
    toappend = {
        "go":go_category,
        "chi2":calculate_chi_square(contingency_table2,expected)
    }
    list_of_go_terms_with_chi2.append(toappend)

for elem in list_of_go_terms_with_chi2:
    print(elem)
