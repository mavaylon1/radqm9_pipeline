import pickle

def read_elements(file: str):
    """
    Read pickle from the file path provided.
    """
    with open(file, 'rb') as io:
        _dict = pickle.load(io)
        return _dict