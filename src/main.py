from src.Coil import Coil
from src.MatchCreator import MatchCreator
import gurobipy as gp
import numpy as np
from gurobipy import GRB


def check_adjacent_matches(match1: tuple, match2: tuple):
    """
            Initializes object coil

            Parameters
            ----------
            :param match1: tuple
                first match in form (coil1: Coil, strand, mode)
            :param match2: tuple
                second match in form (coil2: Coil, strand, mode)

            :returns bool
                True if possible without stringer
    """
    output = True
    coil1 = match1[0]
    mode1 = match1[2]
    coil2 = match2[0]
    mode2 = match2[2]

    ##############################################################
    #            SET PARAMETERS HERE
    ##############################################################
    thickness_threshold = 0.15
    width_threshold = 15
    temperature_threshold = 10
    speed_threshold = 10

    if abs(coil1.thickness - coil2.thickness) > thickness_threshold:
        output = False
    if abs(coil1.width - coil2.width) > width_threshold:
        output = False

    if abs(mode1[0] - mode2[0]) > temperature_threshold:
        output = False
    if abs(mode1[1] - mode2[1]) > speed_threshold:
        output = False

    return output


def compute_A_p_q_k(groups: dict[int, list[tuple[Coil, int, tuple[int, int]]]]) -> list[np.ndarray]:
    """
    Calculate matrix A of compatibility

    :param groups: dict[int, list[tuple[Coil, int, tuple[int, int]]]]

    :return: list[np.ndarray]
        A[k][p, q] = 1 if matches in position p and q in group k are incompatible, 0 otherwise
        p, q real positions from 1... because of dummy initial and end pos
    """
    A = []

    for k in groups.keys():

        group_length = len(groups[k])
        group_matrix = np.zeros((group_length, group_length))
        for p in range(group_length):
            for q in range(group_length):
                if not check_adjacent_matches(groups[k][p], groups[k][q]):
                    group_matrix[p][q] = 1
        A.append(group_matrix)
    return A


def compute_C_k(groups: dict[int, list[tuple[Coil, int, tuple[int, int]]]], dict_of_coils: dict[int, int]) \
        -> dict[int, list[int]]:
    """
    Calculate dict C, each key corresponds to a group, value represents
    the set of coils in group M_k (without repetition)

    :param groups: dict[int, list[tuple[Coil, int, tuple[int, int]]]]
    :param dict_of_coils: dict[int, int]

    :return: dict[int, list[int]]
    dict of pairs group_id -> list of coils in the group
    """
    C = dict()

    for key, value in groups.items():
        coils_without_rep = set()
        for match in value:
            coil_id = match[0].coil_id
            mip_translated_coil_id = list(dict_of_coils.values()).index(coil_id)
            coils_without_rep.add(mip_translated_coil_id)
        C[key] = list(coils_without_rep)
    return C


def compute_D(dict_of_coils: dict[int, int], C: dict[int, list[int]]) -> dict[int, list[int]]:
    """
    Calculate dict D, each key represents coil id (as in input) and value is list of group_ids
    where the coil is contained

    :param dict_of_coils: dict[int, int]
        list of coil objects
    :param C: dict[int, list[int]]
         as computed in the function above (each key corresponds to a group, value represents
    the set of coils in group M_k (without repetition))
    :return: dict[int, list[int]]
    """

    D = dict()
    for c in dict_of_coils.keys():
        D[c] = []
    for key, value in C.items():
        for mip_translated_coil_id in value:
            D[mip_translated_coil_id].append(key)
    return D

def E_i_k(k: int, i: int, groups: dict[int, list[tuple[Coil, int, tuple[int, int]]]], dict_of_coils: dict[int, int])\
                                -> list[int]:
    """
    Computes E_i_k from MIP. E_i_k is a list of positions of matches related to coil i (index as in MIP)
    i group k

    :param k: int
        group id
    :param i: int
        mip translated coil id
    :param groups: dict[int, list[tuple[Coil, int, tuple[int, int]]]]
    :param dict_of_coils: dict[int, int]

    :return: list[int]
    """
    out = []
    for position, match in enumerate(groups[k]):
        coil_id = match[0].coil_id
        for key, value in dict_of_coils.items():
            if value == coil_id:
                mip_translated_coil_id = key
                break
        if mip_translated_coil_id == i:
            out.append(position)
    return out


if __name__ == '__main__':
    list_of_coils: list[Coil] = []


    # reading input file
    file_name: str = 'data_coils.csv'
    with open('../data/' + file_name) as f:
        for count, line in enumerate(f.readlines()):
            # skip first line (header)
            if count == 0:
                continue

            coil_info = line.split(';')

            # Coil(length: int, width: int, thickness: float, min_temp: int, max_temp: int, min_speed: int,
            #                  max_speed: int, coil_id: int )
            list_of_coils.append(Coil(int(float(coil_info[0])), int(coil_info[1]), float(coil_info[2]),
                                      int(float(coil_info[3])), int(float(coil_info[4])), int(float(coil_info[5])),
                                      int(float(coil_info[6])), int(float(coil_info[8]))))

    """
    For MIP we need numeration from 0, this is not always the case with input data.
    list_of_coils has the data from file, dict_of_coils is translation MIP_id -> input coil_id
    """
    dict_of_coils = dict()
    for i, coil in enumerate(list_of_coils):
        dict_of_coils[i] = coil.coil_id
    """
    rules for available strands(lines), a tuple for each line in the form 
                                    (min_allowed_width (cm), max_allowed_width (cm)) 
    """
    strands = ((0, 75), (75, 100), (50, 100))
    mc = MatchCreator(list_of_coils, strands)

    """ 
    forming groups as in Step 2 of Group and greedy pick (sec 5.1)
        1. sort the matches on each line according to some criteria
        2. iterate through the sorted list and check the compatibility of adjacent matches
    """
    first_line_sorted = sorted(mc.get_matches_on_line(1), key=lambda coil: coil[2][1])  # change thickness
    second_line_sorted = sorted(mc.get_matches_on_line(2), key=lambda coil: coil[2][1])  # change thickness
    third_line_sorted = sorted(mc.get_matches_on_line(3), key=lambda coil: coil[2][1])  # change thickness for

    groups_first_line: list[list[tuple[Coil, int, tuple[int, int]]]] = []  # group is an ordered list of matches
    groups_second_line: list[list[tuple[Coil, int, tuple[int, int]]]] = []  # group is an ordered list of matches
    groups_third_line: list[list[tuple[Coil, int, tuple[int, int]]]] = []  # group is an ordered list of matches

    ################################################################################################
    # For usage in MIP all groups has be indexed together as strands (lines) are not included in MIP
    ################################################################################################
    groups: dict[int, list[tuple[Coil, int, tuple[int, int]]]] = dict()  # key index used in MIP, value group

    group_index: int = 0
    group = []
    for i in range(len(first_line_sorted)):
        if i != len(first_line_sorted) - 1:  # exclude the last item to avoid integer out of range error in next step
            if len(group) == 0:  # group has at least one match
                group.append(first_line_sorted[i])
            if check_adjacent_matches(first_line_sorted[i], first_line_sorted[i + 1]):
                group.append(first_line_sorted[i + 1])
            else:
                groups_first_line.append(group)
                groups[group_index] = group
                group_index += 1
                group = [first_line_sorted[i + 1]]
        else:
            groups_first_line.append(group)
            groups[group_index] = group
            last_index_first_line = group_index
            group_index += 1

    group = []  # same now for other lines, indexing continued
    for i in range(len(second_line_sorted)):
        if i != len(second_line_sorted) - 1:  # exclude the last item to avoid integer out of range error in next step
            if len(group) == 0:  # group has at least one match
                group.append(second_line_sorted[i])
            if check_adjacent_matches(second_line_sorted[i], second_line_sorted[i + 1]):
                group.append(second_line_sorted[i + 1])
            else:
                groups_second_line.append(group)
                groups[group_index] = group
                group_index += 1
                group = [second_line_sorted[i + 1]]
        else:
            groups_second_line.append(group)
            groups[group_index] = group
            last_index_second_line = group_index
            group_index += 1

    group = []

    for i in range(len(third_line_sorted)):
        if i != len(third_line_sorted) - 1:  # exclude the last item to avoid integer out of range error in next step
            if len(group) == 0:  # group has at least one match
                group.append(third_line_sorted[i])
            if check_adjacent_matches(third_line_sorted[i], third_line_sorted[i + 1]):
                group.append(third_line_sorted[i + 1])
            else:
                groups_third_line.append(group)
                groups[group_index] = group
                group_index += 1
                group = [third_line_sorted[i + 1]]
        else:
            groups_second_line.append(group)
            groups[group_index] = group
            last_index_third_line = group_index

    """ Compute parameters needed for MIP """
    A = compute_A_p_q_k(groups)
    C = compute_C_k(groups, dict_of_coils)
    D = compute_D(dict_of_coils, C)

    K = len(A)  # number of groups
    N = len(list(D))    # number of coils





    m = gp.Model("GOPH")


    Y = []
    for k, group in enumerate(A):
        n = A[k].shape[0]
        Y.append(m.addVars(n + 1, n + 1, vtype=GRB.BINARY, name="y"))   # one starting and one ending dummy position

    Z = m.addVars(K, vtype=GRB.BINARY, name="z")

    X = m.addVars(N, K, vtype=GRB.BINARY, name="x")

    m.setObjective(sum(sum(sum(A[k][p, q] * Y[k][p, q] for q in range(len(groups[k])) if q > p)
                           for p in range(len(groups[k]))) for k in range(K)) +
                   sum(Z[k] for k in range(K)))

    m.addConstrs(sum(Y[k][p, q] for q in range(1, len(groups[k]) - 1) if q > p) ==
                 sum(Y[k][p, q] for q in range(1, len(groups[k]) - 1) if q < p) for k in range(K)
                 for p in range(len(groups[k])))

    m.addConstrs(sum(Y[k][0, q] for q in range(1, len(groups[k]) + 1)) == 1 for k in range(K))

    m.addConstrs(sum(Y[k][p, len(groups[k])] for p in range(len(groups[k]))) == 1
                 for k in range(K))
    m.addConstrs(X[i, k] == sum(sum(Y[k][p, q] for q in range(len(groups[k])) if q > p)
                                    for p in E_i_k(k, i, groups, dict_of_coils)) for i in C[k]
                                    for k in range(K))

    m.addConstrs(X[i, k] <= Z[k] for k in range(K) for i in C[k])

    m.addConstrs(sum(X[i, k] for k in D[i]) == 1 for i in range(N))
    m.optimize()

    result = m.getObjective().getValue()

