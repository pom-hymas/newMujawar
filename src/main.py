

# Press the green button in the gutter to run the script.
from src.Coil import Coil
from src.MatchCreator import MatchCreator

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
        print("thickness")
        output = False
    if abs(coil1.width - coil2.width) > width_threshold:
        output = False

    if abs(mode1[0] - mode2[0]) > temperature_threshold:
        output = False
    if abs(mode1[1] - mode2[1]) > speed_threshold:
        output = False

    return output

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
    first_line_sorted = sorted(mc.get_matches_on_line(1), key=lambda coil: coil[0].thickness)  # change thickness
    second_line_sorted = sorted(mc.get_matches_on_line(2), key=lambda coil: coil[0].thickness)  # change thickness
    third_line_sorted = sorted(mc.get_matches_on_line(3), key=lambda coil: coil[0].thickness)  # change thickness for

    groups_first_line: list[list[tuple[Coil, int, tuple[int, int]]]] = []  # group is an ordered list of matches
    groups_second_line: list[list[tuple[Coil, int, tuple[int, int]]]] = []  # group is an ordered list of matches
    groups_third_line: list[list[tuple[Coil, int, tuple[int, int]]]] = []  # group is an ordered list of matches

    for i in range(len(first_line_sorted)):
        if i != len(first_line_sorted) - 1:
            pass
            #print(check_adjacent_matches(first_line_sorted[i], first_line_sorted[i+1]))




