# -*- coding: utf-8 -*-

def read(ws, filename):
    """ Reads an online Hitran catalog of the latest ARTS specifications

    These specifications are:

    Field separator: [tab], Line endings: Windows (CR LF)

    .par_line
    qns'
    qns''

    Returns are reference to the generated abs_lines
    
    Parameters:
        ws (Workspace): A pyarts workspace
        filename (str): The absolute or relative path to a Hitran file
    
    Return:
        The read abs_lines as a pyarts workspace variable
    """
    ws.ReadHITRAN(filename=filename,
      normalization_option="SFS",
      hitran_type="Online")
    ws.abs_linesSetEmptyBroadeningParametersToEmpty()

    return ws.abs_lines
