"""
ezdxf_test.py
prints out the points-list of the first LWPolyline entity
"""
import sys
import ezdxf

try:
    doc = ezdxf.readfile(str(sys.argv[1]))
except IOError:
    print(f'Not a DXF file or a generic I/O error.')
    sys.exit(1)
except ezdxf.DXFStructureError:
    print(f'Invalid or corrupted DXF file.')
    sys.exit(2)

msp = doc.modelspace()
# take first LWPolyline, 'first' was introduced with v0.10
line = msp.query('LWPOLYLINE').first
with line.points('xyseb') as points:
    p_str = []
    for point in points:
        p_str.append(tuple(float(v) for v in point))
    print(p_str)