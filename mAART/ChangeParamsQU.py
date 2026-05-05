import argparse
import re
import sys

def update_param(content, param_name, new_value, is_int=False):
    """Updates a parameter value. If is_int is True, forces integer format."""
    # Convert to float first to handle cases like "20.0" being passed as a string
    val = float(new_value)
    if is_int:
        val_str = str(int(val))
    else:
        # Use g to keep 0.02 as 0.02 but allow float precision
        val_str = f"{val:g}"
        if "." not in val_str and not is_int:
            val_str = f"{val:.2f}" # Keep spin as 0.02 instead of 0.020000

    pattern = r"^\s*" + re.escape(param_name) + r"\s*=\s*[\d\.\-eE]+"
    replacement = f"{param_name}={val_str}"

    return re.sub(pattern, replacement, content, flags=re.MULTILINE)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--a', type=str)
    parser.add_argument('--i', type=str)
    parser.add_argument('--Br', type=str)
    parser.add_argument('--Bth', type=str)
    parser.add_argument('--Bphi', type=str)
    parser.add_argument('--sub_kep', type=str)
    parser.add_argument('--betar', type=str)
    parser.add_argument('--betaphi', type=str)

    args = parser.parse_args()

    with open('params.py', 'r') as f:
        content = f.read()

    if args.a:
        content = update_param(content, 'spin_case', args.a)
    if args.i:
        # FORCE INTEGER HERE
        content = update_param(content, 'i_case', args.i, is_int=True)
    if args.Br:
        content = update_param(content, 'Br', args.Br)
    if args.Bth:
        content = update_param(content, 'Bth', args.Bth)
    if args.Bphi:
        content = update_param(content, 'Bphi', args.Bphi)
    if args.sub_kep:
        content = update_param(content, 'sub_kep', args.sub_kep)
    if args.betar:
        content = update_param(content, 'betar', args.betar)
    if args.betaphi:
        content = update_param(content, 'betaphi', args.betaphi)

    with open('params.py', 'w') as f:
        f.write(content)

if __name__ == "__main__":
    main()
