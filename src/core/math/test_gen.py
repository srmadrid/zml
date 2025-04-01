# Temporary script to generate Zig test cases from a test vector file
import sys
import os
from collections import defaultdict

def main():
    # Make it work with multiple test vector files?
    if len(sys.argv) < 3:
        print('Usage: python generate_tests.py <test_vector_file> <output_file>')
        return

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    with open(input_file, 'r') as f:
        lines = f.readlines()

    # {function_name: {'is_complex': bool, 'tests': {type: [(args, expected)]}}}
    test_cases = defaultdict(lambda: {
        'is_complex': False,
        'tests': defaultdict(list)
    })

    for line in lines:
        line = line.strip()
        if not line.startswith('= '):
            continue
            
        parts = line.split()
        if len(parts) < 9:
            continue

        func_name = parts[1]
        rounding = parts[2]
        fmt = parts[3]
        args = parts[4]
        expected = parts[6]

        if rounding != 'tonearest':
            continue

        # Determine function type based on name prefix
        is_complex = func_name.startswith('c')
        test_cases[func_name]['is_complex'] = is_complex

        # Get type mapping
        type_map = {
            'binary32': 'f32',
            'binary64': 'f64',
            'intel96': 'f80',
            'binary128': 'f128'
        }
        if fmt not in type_map:
            continue
            
        zig_type = type_map[fmt]
        
        if args == 'plus_infty':
            args = f'std.math.inf({zig_type})'
        elif args == 'minus_infty':
            args = f'-std.math.inf({zig_type})'
        
        if expected == 'plus_infty':
            expected = f'std.math.inf({zig_type})'
        elif expected == 'minus_infty':
            expected = f'-std.math.inf({zig_type})'

        # Store test case with proper formatting
        if is_complex:
            # Complex function: single argument with real/imag parts
            if len(args) != 2:
                continue
            test_args = f'types.c{zig_type}.init({args[0]}, {args[1]})'
        else:
            # Real function: multiple arguments
            if args == f'std.math.inf({zig_type})' or args == f'-std.math.inf({zig_type})':
                test_args = args
            else:
                test_args = f'@as({zig_type}, {args})'

        test_cases[func_name]['tests'][zig_type].append((test_args, expected))

    # Generate Zig test code
    zig_code = []

    for func_name, data in test_cases.items():
        zig_code.append(f'test {func_name} {{')
        
        for zig_type in ['f32', 'f64', 'f80', 'f128']:
            tests = data['tests'].get(zig_type, [])
            if not tests:
                continue
                
            # zig_code.append(f'    // {zig_type}')
            for args, expected in tests:
                test_line = f'    try std.testing.expectEqual('
                test_line += f'{expected}, {func_name}({args}));'
                zig_code.append(test_line)
            
            zig_code.append('')

        zig_code.append('}\n')

    with open(output_file, 'w') as f:
        f.write('\n'.join(zig_code))

    print(f'Generated {output_file} with:')
    for func_name, data in test_cases.items():
        print(f' - {func_name}: {sum(len(tests) for tests in data['tests'].values())} tests')

if __name__ == '__main__':
    main()