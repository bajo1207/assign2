def get_list_elem(list, index):
    try:
        return list[index].rstrip()
    except IndexError:
        return None

def create_list(filename):
    lst = []
    file = open(filename, "r")
    for line in file:
        lst.append(line.rstrip())
    file.close()
    return lst

def calculate_chi_expected_values(contingency_table):
    for n in contingency_table:
        n.append(sum(n))
    col_totals = [ sum(x) for x in zip(*contingency_table) ]
    contingency_table.append(col_totals)
    expected = []
    for r in range(0, len(contingency_table) - 1):
        row = []
        for c in range(0, len(contingency_table) - 1):
            row.append((contingency_table[r][-1] * contingency_table[-1][c])/contingency_table[-1][-1])
        expected.append(row)
    return expected

def calculate_chi_square(observed, expected):
    total = 0
    for (a,b) in zip(observed, expected):
        for (o,e) in zip(a,b):
            total += (o-e)**2/e
    return total

            
