import yaml
from os.path import join, abspath, dirname
from pprint import pprint

supplement_dir = dirname(abspath(__file__))

with open(join(supplement_dir, 'scanpy.yml'), 'r') as f:
    sc_pl_data = yaml.safe_load(f)

with open(join(supplement_dir, 'spatialdata-plot.yml'), 'r') as f:
    sdata_pl_data = yaml.safe_load(f)

#pprint(sc_pl_data)
#pprint(sdata_pl_data)

output = ""

PARAM_STATES = [
    'supported',
    'partially supported',
    'not currently supported',
    'lacks description',
    'N/A',
]

output += "# Scanpy plotting functions\n\n"

for func_name, func_data in sc_pl_data.items():
    output += f"## sc.pl.{func_name}\n\n"
    
    func_info = func_data['function']
    is_implemented = func_info.get('implemented')
    if not is_implemented:
        output += "_Not yet implemented._\n\n"
        continue

    param_info = func_data['parameters']
    params_by_state = {state: [] for state in PARAM_STATES}
    for param_name, param_state in param_info.items():
        params_by_state[param_state].append(param_name)

    for param_state, params in params_by_state.items():
        num_params = len(params)
        total_params = sum(len(p) for p in params_by_state.values())
        pct_params = (num_params / total_params * 100) if total_params > 0 else 0
        if len(params) > 0:
            output += f"### Parameter state: {param_state}\n\n"
            for param in params:
                output += f"- `{param}`\n"
            output += "\n"

output += "# SpatialData-Plot plotting functions\n\n"

for func_name, func_data in sdata_pl_data.items():
    output += f"## sdata.pl.{func_name}\n\n"
    
    func_info = func_data['function']
    is_implemented = func_info.get('implemented')
    if not is_implemented:
        output += "_Not yet implemented._\n\n"
        continue

    param_info = func_data['parameters']
    params_by_state = {state: [] for state in PARAM_STATES}
    for param_name, param_state in param_info.items():
        params_by_state[param_state].append(param_name)

    for param_state, params in params_by_state.items():
        num_params = len(params)
        total_params = sum(len(p) for p in params_by_state.values())
        pct_params = (num_params / total_params * 100) if total_params > 0 else 0
        if len(params) > 0:
            output += f"### Parameter state: {param_state}\n\n"
            for param in params:
                output += f"- `{param}`\n"
            output += "\n"

print(output)
