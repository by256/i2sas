import numpy as np


def parse_datathief(data_path):
    """
    Parse .txt file produced by datathief. Returns np.array
    """
    data = []
    with open(data_path, 'r') as f:
        for line in f.readlines()[1:]:
            line = line.strip('\n')
            line = line.split(', ')
            line = [float(x) for x in line]
            data.append(line)
    data = np.array(data)
    sort_idx = np.argsort(data[:, 0])
    data = data[sort_idx]
    return data
