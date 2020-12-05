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
    ct = []
    for row in contingency_table:
        ct.append(row.copy())
    for n in ct:
        n.append(sum(n))
    col_totals = [ sum(x) for x in zip(*ct) ]
    ct.append(col_totals)
    expected = []
    for r in range(0, len(ct) - 1):
        row = []
        for c in range(0, len(ct) - 1):
            row.append((ct[r][-1] * ct[-1][c])/ct[-1][-1])
        expected.append(row)
    return expected

def calculate_chi_square(observed, expected):
    total = 0
    for (o,e) in zip(observed, expected):
        for (o,e) in zip(o,e):
            total += (o-e)**2/e
    return total

            
