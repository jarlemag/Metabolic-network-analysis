#xmltodict.py


def recursive_dict(element):
    return element.tag, \
        dict(map(recursive_dict, element)) or element.text
