import json

def save_json(data, path, indent=2):
    """
    Save a dictionary as a pretty-printed JSON file.
    """
    with open(path, 'w') as f:
        json.dump(data, f, indent=indent)


def pretty_print_json(data, indent=2):
    """
    Return a pretty-printed JSON string.
    """
    return json.dumps(data, indent=indent) 