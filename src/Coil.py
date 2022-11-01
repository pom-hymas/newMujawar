class Coil:
    """
    A class representing coils

    Attributes
    ----------
    coil_id: int
    length: int
        length in m
    width: int
        width in cm
    thickness: float
        thickness in cm
    min_temp: int
        minimal allowed temp for the coil
    max_temp: int
        maximal allowed temp for the coil
    min_speed: int
        min strip speed allowed for the coil
    max_speed: int
        max strip speed allowed for the coil
    """
    def __init__(self, length: int, width: int, thickness: float, min_temp: int, max_temp: int, min_speed: int,
                 max_speed: int, coil_id: int):
        """
        Initializes object coil

        Parameters
        ----------
        :param length: int
            length in
        :param width: int
            width in cm
        :param thickness: float
            thickness in cm
        :param coil_id: int
            coil id
        """
        self.coil_id = coil_id
        self.length = length
        self.width = width
        self.thickness = thickness
        self.min_temp = min_temp
        self.max_temp = max_temp
        self.min_speed = min_speed
        self.max_speed = max_speed

        temp_modes = [temp_mode for temp_mode in range(min_temp, max_temp + 1, 10)]
        speed_modes = [speed_mode for speed_mode in range(min_speed, max_speed + 1, 10)]

        self.modes = [(temp, speed) for temp in temp_modes for speed in speed_modes]

    def __str__(self):
        out = "Coil ID: {}," \
              "length: {}" \
              "width: {}" \
              "thickness: {}".format(self.coil_id, self.length, self.width, self.thickness)
        return out