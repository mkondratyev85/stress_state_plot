from plane import Plane

def parse_line(line):
    line = line.strip()
    line = line.split(" ")
    direction = float(line[0])
    dip = float(line[1])
    return direction, dip

def load_fractures(filename):
    fractures = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line_number, line in enumerate(lines):
            try: 
                dr, dp = parse_line(line)
                fractures.append(Plane(dr, dp))
            except ValueError as e:
                print("Error while reading line number {line_number}")
                print(e)
    return fractures
