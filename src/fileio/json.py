import json
from functools import partial

def json_parse(fileobj, decoder=json.JSONDecoder(), buffersize=2048):
    """Function taken from:
    https://stackoverflow.com/a/21709058


    Arguments:
        fileobj {[type]} -- [description]

    Keyword Arguments:
        decoder {[type]} -- [description] (default: {JSONDecoder()})
        buffersize {int} -- [description] (default: {2048})
    """
    buffer = ''
    for chunk in iter(partial(fileobj.read, buffersize), ''):
         buffer += chunk
         while buffer:
             try:
                 result, index = decoder.raw_decode(buffer)
                 yield result
                 buffer = buffer[index:].lstrip()
             except ValueError:
                 # Not enough data to decode, read more
                 break


def read_json(fname):
    """Read the VMS json

    Arguments:
        fname {str} -- filename
    """
    with open(fname, 'r') as fopen:
        for data in json_parse(fopen):
            # process object
            res = data

    return res
