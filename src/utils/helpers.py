import random
from typing import List

def readlistnum(x, nmax=None):
    result = []
    for part in x.split(','):
        if '-' in part:
            a, b = part.split('-')
            if a == "": 
                a = 1
            else:       
                a = int(a)
            if b == "":
                if nmax:
                    b = int(nmax)
                else:
                    raise NotImplementedError()
            elif int(b) > nmax:
                raise IndexError
            else:
                b = int(b)
            result.extend(range(a, b + 1))
        else:
            if int(part) > nmax:
                raise IndexError
            else:
                a = int(part)
                result.append(a)
    return result

def random_colors(ncolors: int) -> List[str]:
    """Return a list of n random colors.

    Args:
        ncolors: Number of colors to generate
        
    Returns:
        List of hex color strings
    """
    color = ["#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(ncolors)]
    return color
