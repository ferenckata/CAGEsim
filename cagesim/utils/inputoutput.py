'''Input/output class for file writing and reading methods'''
# TODO: logging.log
import time


class InputOutput:
    """Class with basic file reading / writing functionalities"""

    @classmethod
    def list_to_file(cls, data_list: list, filename: str) -> None:
        '''Export data in list format to a txt file
        One entry = one line
        Parameters
        ----------
        data_list:
            any type of data in a list format
        filename:
            output file name (without extension)
        '''
        timestr = time.strftime("%Y%m%d-%H%M%S")
        with open(filename + '_' + timestr + '.txt', 'w', encoding='utf-8') as out_file:
            for line in data_list:
                out_file.write(f"{line}\n")
