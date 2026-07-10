"""
Generic configuration parsing and validation logic.
"""
import sys
import yaml

def get_provided_cli_flags():
    """
    Extract a simple list of explicitly provided CLI flags.
    TODO: There may be a more robust approach than scanning sys.argv directly, 
    such as hooking into argparse.Action or using parser.parse_known_args.
    """
    provided = set()
    for arg in sys.argv:
        if arg.startswith('-'):
            # Handle arguments like --max-iter=2
            flag = arg.split('=')[0]
            provided.add(flag)
    return provided

def merge_config(args, config_mapping):
    """
    Read the YAML config (if args.config is provided), validate keys and types,
    check for mutual exclusion with explicitly provided CLI arguments, 
    and merge the config values into the args namespace.
    
    Returns True if successful, False if an error occurred.
    """
    if not getattr(args, 'config', None):
        return True
        
    try:
        with open(args.config, 'r') as f:
            config_data = yaml.safe_load(f)
    except Exception as e:
        print(f"ERROR: Could not read config file {args.config}: {e}", file=sys.stderr)
        return False
        
    if config_data is None:
        config_data = {}
        
    # Check for unknown keys
    unknown_keys = set(config_data.keys()) - set(config_mapping.keys())
    if unknown_keys:
        first_unknown = list(unknown_keys)[0]
        # Match error message from tests: ERROR: Invalid refine config 'unknown.yaml': unexpected config option 'unknown'
        command_name = sys.argv[1] if len(sys.argv) > 1 else 'command'
        print(f"ERROR: Invalid {command_name} config '{args.config}': unexpected config option '{first_unknown}'", file=sys.stderr)
        return False

    provided_cli_flags = get_provided_cli_flags()
    conflicts = []

    for key, value in config_data.items():
        mapping = config_mapping[key]
        expected_type = mapping.get('type')
        
        # Type checking
        # TODO: jsonschema may be more robust for validation in the future.
        if expected_type and not isinstance(value, expected_type):
            print(f"ERROR: Invalid type for '{key}' in config file '{args.config}'. Expected {expected_type}, got {type(value)}", file=sys.stderr)
            return False
            
        # Mutual exclusion
        cli_flags = mapping.get('cli', [])
        conflicting_flags = [flag for flag in cli_flags if flag in provided_cli_flags]
        if conflicting_flags:
            conflicts.extend(conflicting_flags)
        else:
            setattr(args, key, value)
            
    if conflicts:
        print(f"ERROR: --config cannot be used with these CLI options: {', '.join(conflicts)}", file=sys.stderr)
        return False
        
    return True
