'''' Practical Computing for Bioinformatics: Assignment 2

MAIN TASK: Calculate the enrichment of specific GO terms in a list of selected genes and output the top 5 results

authors: Joachim Bache-Mathiesen, Antonietta Salerno
date: 8/11/2020

'''

from GOenrich_utils import *

# Read GO annotations and create a list of dictionaries, where every element is a line
gotext = open("GOannotations.txt", "r")
goannotations = [] # list of dictionaries with the split lines
allgenes = [] # all genes having a GO term (coding genes)
unique_set_of_goterms = set() # unique set of all GO IDs
next(gotext) # skip first line: columns titles
for line in gotext:
    content_array = line.split(maxsplit=2)
    item = {
        "gene_stab_id":content_array[0],
        "go_term":content_array[1],
        "term_name":get_list_elem(content_array, 2)
    } 
    goannotations.append(item)
    
    ''' DISCLAIM:
     We decided to directly take under examination those genes that for sure have a GO term associated.
     to find the genes' enrichment for a specific GO term in order to mitigate computational efforts.
     This choice is also driven by the fact that we noticed that there are gene IDs that exist in the GOannotations file but are not contained in allGenes, 
     this does not only provide a less reliable statistic but also implies the computational problem of the division by 0.
     '''
    
    allgenes.append(item["gene_stab_id"]) # Store all genes found in the GOannotation file
    unique_set_of_goterms.add(item['go_term']) 
unique_list_of_goterms = list(unique_set_of_goterms) # Create a list containing all the GO terms
gotext.close()

selected_genes = create_list("IDS2.txt") # Read the file containing the selected genes
allgenes_boolean_selected = {}
for gene in allgenes:
    sel = gene in selected_genes
    allgenes_boolean_selected[gene] = sel # We create a dictionary having as keys all genes and as value whether they are found in the selection or not


"""
In order to see if there are more occurrences of a GO term in the gene selection than expected 
we need to look at how many genes are associated to that category and at the same time find out if those gene are present in the selection or not.

We generate a new file containing for each GO term, its ID, name and chi-squared value.
This last is an indicator of how much that GO term is overrepresented in the selection.
"""

f = open('go_with_chi2.txt', 'w') 
list_of_go_terms_with_chi2 = []
for go_term in unique_list_of_goterms:
    go_term_genes = [] # list of genes associated to a specific GO term
    contingency_table = [[0,0],[0,0]]
    #Find all gene IDs that belong to a given GO term
    for go in goannotations:
        if go['go_term'] == go_term:
            go_term_genes.append(go['gene_stab_id'])
            go_term_name = go["term_name"] # report their names
            
    # We have to look at all genes to fill out a contingency table for each GO term: needed to compute the chi-squared
    for gene, sel in allgenes_boolean_selected.items():
        go = gene in go_term_genes
        if go and sel:
            contingency_table[0][0] += 1
        if go and (not sel):
            contingency_table[0][1] += 1
        if (not go) and sel:
            contingency_table[1][0] += 1
        if (not go) and (not sel): 
            contingency_table[1][1] += 1
    expected = calculate_chi_expected_values(contingency_table)
    chi2 = calculate_chi_square(contingency_table,expected)
    f.write(f"{go_term}    {go_term_name}     {chi2}\n")
    toappend = {
        "go_term":go_term,
        "go_term_name": go_term_name,
        "chi2":chi2
    }
    list_of_go_terms_with_chi2.append(toappend)
 
# We output the 5 most overrepresented GO terms in the gene selection together with their chi-squared statistic
newlist = sorted(list_of_go_terms_with_chi2, key=lambda k: k['chi2'])
for item in newlist[-5:]:
    print(item)

