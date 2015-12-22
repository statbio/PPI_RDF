import sys

def main():
    f = open(sys.argv[1])
    lines = f.readlines()
    f.close()

    edited_lines = []
    for i, line in enumerate(lines):
        if len(line) < 2:
            continue
        if line[0] is '>':
            new_line = lines[i-1].strip('\n') + line
            edited_lines.append(new_line)
        elif line[-2] is '.':
            edited_lines.append(line)


    output = open('BIOGRID_fixed.nt', 'w')
    for nl in edited_lines:
        output.write(nl)
    output.close()


if __name__ == "__main__":
    main()
