"""Test that all plot functions follow the signature convention.

All plot functions should have:
1. First parameter named 'data' (positional)
2. All other parameters are keyword-only (after *)
"""

import ast
import inspect
from pathlib import Path
import importlib.util


def get_plot_modules():
    """Get all plot module files from pyarts3.plots."""
    plots_dir = Path(__file__).parent.parent.parent / "src" / "pyarts3" / "plots"
    
    # Get all .py files except __init__.py
    plot_files = [f for f in plots_dir.glob("*.py") 
                  if f.name != "__init__.py" and not f.name.startswith("test")]
    
    return plot_files


def check_function_signature(func_node, filename):
    """Check if a function has the correct signature.
    
    Returns tuple: (is_valid, error_message)
    """
    args = func_node.args
    
    # Get positional arguments (before *)
    positional_args = args.args
    
    # Get keyword-only arguments (after *)
    kwonly_args = args.kwonlyargs
    
    # Check if there are exactly 1 positional argument
    if len(positional_args) != 1:
        return False, f"{filename}:{func_node.lineno} - Function '{func_node.name}' has {len(positional_args)} positional args, expected 1"
    
    # Check if the first positional argument is named 'data'
    first_arg = positional_args[0]
    if first_arg.arg != 'data':
        return False, f"{filename}:{func_node.lineno} - Function '{func_node.name}' first parameter is '{first_arg.arg}', expected 'data'"
    
    # All checks passed
    return True, None


def test_plot_function_signatures():
    """Test that all plot functions have correct signatures."""
    errors = []
    
    plot_files = get_plot_modules()
    
    for plot_file in plot_files:
        # Parse the file
        with open(plot_file, 'r') as f:
            try:
                tree = ast.parse(f.read(), filename=str(plot_file))
            except SyntaxError as e:
                errors.append(f"Syntax error in {plot_file.name}: {e}")
                continue
        
        # Find all function definitions
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef):
                # Only check functions starting with 'plot'
                if not node.name.startswith('plot'):
                    continue
                
                # Check the function signature
                is_valid, error_msg = check_function_signature(node, plot_file.name)
                if not is_valid:
                    errors.append(error_msg)
    
    # Report errors
    if errors:
        error_report = "\n".join(errors)
        raise AssertionError(f"Found {len(errors)} signature violations:\n{error_report}")


def test_plot_function_signatures_detailed():
    """Detailed test with examples of correct vs incorrect signatures."""
    plot_files = get_plot_modules()
    
    results = {
        'correct': [],
        'incorrect': []
    }
    
    for plot_file in plot_files:
        with open(plot_file, 'r') as f:
            try:
                tree = ast.parse(f.read(), filename=str(plot_file))
            except SyntaxError:
                continue
        
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef):
                # Only check functions starting with 'plot'
                if not node.name.startswith('plot'):
                    continue
                
                is_valid, error_msg = check_function_signature(node, plot_file.name)
                
                # Build signature string for display
                args = node.args
                sig_parts = []
                
                # Positional args
                for arg in args.args:
                    sig_parts.append(arg.arg)
                
                # Add * separator if there are kwonly args
                if args.kwonlyargs:
                    sig_parts.append('*')
                    for arg in args.kwonlyargs:
                        default_marker = '=...' if any(d for d in node.args.kw_defaults) else ''
                        sig_parts.append(f"{arg.arg}{default_marker}")
                
                sig_str = f"{node.name}({', '.join(sig_parts)})"
                
                if is_valid:
                    results['correct'].append((plot_file.name, sig_str))
                else:
                    results['incorrect'].append((plot_file.name, sig_str, error_msg))
    
    # Print summary
    print(f"\n{'='*70}")
    print(f"Plot Function Signature Test Summary")
    print(f"{'='*70}")
    print(f"Correct: {len(results['correct'])}")
    print(f"Incorrect: {len(results['incorrect'])}")
    
    if results['incorrect']:
        print(f"\n{'='*70}")
        print("INCORRECT SIGNATURES:")
        print(f"{'='*70}")
        for filename, sig, error in results['incorrect']:
            print(f"\n{filename}")
            print(f"  {sig}")
            print(f"  Error: {error}")
    
    if results['correct']:
        print(f"\n{'='*70}")
        print("CORRECT SIGNATURES (sample):")
        print(f"{'='*70}")
        for filename, sig in results['correct'][:5]:
            print(f"{filename}: {sig}")
        if len(results['correct']) > 5:
            print(f"... and {len(results['correct']) - 5} more")
    
    # Fail if there are incorrect signatures
    assert len(results['incorrect']) == 0, f"Found {len(results['incorrect'])} functions with incorrect signatures"


if __name__ == "__main__":
    print("Running plot function signature tests...")
    try:
        test_plot_function_signatures()
        print("All plot functions have correct signatures!")
    except AssertionError as e:
        print(f"Test failed:\n{e}")
        exit(1)
    
    # Run detailed test for summary
    print("\nRunning detailed signature analysis...")
    try:
        test_plot_function_signatures_detailed()
    except AssertionError as e:
        print(f"\nDetailed test failed:\n{e}")
        exit(1)
