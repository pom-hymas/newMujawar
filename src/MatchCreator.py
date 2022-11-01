from src.Coil import Coil


class MatchCreator:
    """
        A class representing matches according to Mujawar et al. 2012


        Attributes
        ----------
        c: Coil
        current coil for wich matches are created
        number_of_lines: int
            total number of available lines
        suitable_lines: list[int]
            list of lines suitable for the coil c

        """

    def __init__(self, coils_list: list[Coil], strands: tuple[tuple[int, int]]):
        """
                Initializes object coil

                Parameters
                ----------
                :param c: Coil
                    Coil that has to be translated into match
                :param strands: tuple[tuple[int, int]]
                    tuple with each element characterizing a corresponding strand
                    (also as tuple in form (minimal allowed width, maximal allowed width))
        """
        matches = []

        # find out which lines are suitable for the coil

        for line_nr, widths in enumerate(strands):
            for c in coils_list:
                if widths[0] <= c.width < widths[1]:
                    for mode in c.modes:

                        matches.append((c, line_nr + 1, mode))

        self.matches = matches

    def get_matches(self):
        return self.matches

    def get_matches_on_line(self, n: int):
        return [match for match in self.matches if match[1] == n]

    def eliminate_matches_corresponding_to_coil(self, coil_id):
        new_matches = [match for match in self.matches if match[0].coil_id != coil_id]
        self.matches = new_matches
