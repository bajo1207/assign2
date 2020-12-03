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
allgenes = []
for item in goannotations:
    allgenes.append(item["gene_stab_id"])

#Read ids2
ids2genes = create_list("IDS2.txt")

#Here we create a list of all the go terms
unique_set_of_goterms = set()
for goterm in goannotations:
    unique_set_of_goterms.add(goterm['go_term'])
unique_list_of_goterms = list(unique_set_of_goterms)
print(len(unique_list_of_goterms))

allgenes_boolean_selected = {}
for gene in allgenes:
    sel = gene in ids2genes
    allgenes_boolean_selected[gene] = sel

"""
In order to see if there are more occurances of a GO category in the gene selection 
than expected we need to look at how many genes and at the same time find out if this gene was selected or not
"""
f = open('go_with_chi2.csv', 'w')
list_of_go_terms_with_chi2 = []
for go_category in unique_list_of_goterms:
    go_ids = []
    contingency_table2 = [[0,0],[0,0]]
    #Find all gene ids that belong to a given go_category
    for go in goannotations:
        if go['go_term'] == go_category:
            go_ids.append(go['gene_stab_id'])
    # We have to look at all genes to compare gene to selected to fill out contingency table
    for gene, sel in allgenes_boolean_selected.items():
        go = gene in go_ids
        if go and sel:
            contingency_table2[0][0] += 1
        if go and (not sel):
            contingency_table2[0][1] += 1
        if (not go) and sel:
            contingency_table2[1][0] += 1
        if (not go) and (not sel): 
            contingency_table2[1][1] += 1
    expected = calculate_chi_expected_values(contingency_table2)
    chi2 = calculate_chi_square(contingency_table2,expected)
    f.write(f"{go_category}    {chi2}\n")
    toappend = {
        "go":go_category,
        "chi2":chi2
    }
    list_of_go_terms_with_chi2.append(toappend)

newlist = sorted(list_of_go_terms_with_chi2, key=lambda k: k['chi2'])
for item in newlist[-5:]:
    print(item)
