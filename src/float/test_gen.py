# Temporary script to generate Zig test cases from a test vector file
import sys
import os
from collections import defaultdict

def main():
    # Make it work with multiple test vector files?
    if len(sys.argv) < 4:
        print('Usage: python generate_tests.py <arg_count> <test_vector_file> <output_file>')
        return

    arg_count = int(sys.argv[1])
    input_file = sys.argv[2]
    output_file = sys.argv[3]

    with open(input_file, 'r') as f:
        lines = f.readlines()

    # {function_name: {'is_complex': bool, 'tests': {type: [(args, expected)]}}}
    test_cases = defaultdict(lambda: {
        'is_complex': False,
        'tests': defaultdict(list)
    })
    
    is_complex = lines[0].startswith('c') and lines[0][:3] != 'cos' and lines[0][:4] != 'cbrt'
    
    if is_complex:
        arg_count *= 2

    for line in lines:
        line = line.strip()
        if not line.startswith('= '):
            continue
            
        parts = line.split()
        if len(parts) < 7:
            continue

        if is_complex:
            func_name = parts[1][1:]
        else:
            func_name = parts[1]
        
        rounding = parts[2]
        fmt = parts[3]

        if rounding != 'tonearest':
            continue

        test_cases[func_name]['is_complex'] = is_complex
        
        args = parts[4 : 4 + arg_count]
        
        is_complex_output = parts[5 + arg_count + 1][:2] in ['0x', '-0', 'pl', 'mi']
        
        if is_complex_output:
            expected = parts[5 + arg_count: 5 + arg_count + 2]
        else:
            expected = parts[5 + arg_count]

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
        
        for i in range(arg_count):
            if args[i] == 'plus_infty':
                args[i] = f'std.math.inf({zig_type})'
            elif args[i] == 'minus_infty':
                args[i] = f'-std.math.inf({zig_type})'
        
        if is_complex_output:
            for i in range(2):
                if expected[i] == 'plus_infty':
                    expected[i] = f'std.math.inf({zig_type})'
                elif expected[i] == 'minus_infty':
                    expected[i] = f'-std.math.inf({zig_type})'
        else:
            if expected == 'plus_infty':
                expected = f'std.math.inf({zig_type})'
            elif expected == 'minus_infty':
                expected = f'-std.math.inf({zig_type})'
                
        if is_complex_output:
            expected = f'c{zig_type}.init({expected[0]}, {expected[1]})'

        # Store test case with proper formatting
        if is_complex:
            # Complex function:
            if len(args) % 2 != 0:
                continue
            test_args = ''
            for i in range(0, len(args), 2):
                if i > 0:
                    test_args += ', '
                test_args += f'c{zig_type}.init({args[i]}, {args[i + 1]})'
        else:
            # Real function: multiple arguments
            test_args = ''
            for i in range(len(args)):
                if i > 0:
                    test_args += ', '
                if args[i] == f'std.math.inf({zig_type})' or args[i] == f'-std.math.inf({zig_type})':
                    test_args += args[i]
                else:
                    test_args += f'@as({zig_type}, {args[i]})'

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