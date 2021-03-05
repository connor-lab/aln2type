from yaml import safe_load

def read_scheme_yaml(scheme_yaml_path):
    with open(scheme_yaml_path, 'r') as y:
        return safe_load(y)
