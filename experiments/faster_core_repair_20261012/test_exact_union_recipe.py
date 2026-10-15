"""Actual CGAL set semantics, not an optimizer/admission test.

Catches substituting difference for union and omitting the fixed outer shell.
"""
from fractions import Fraction as Q
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

FACES = [[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],
         [1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]]

def box(lo, hi):
    x,y,z=lo; X,Y,Z=hi
    return dict(vertices=[[str(Q(v).numerator)+'/'+str(Q(v).denominator) for v in p]
        for p in [(x,y,z),(X,y,z),(X,Y,z),(x,Y,z),(x,y,Z),(X,y,Z),(X,Y,Z),(x,Y,Z)]],faces=FACES)

def exact_volume(packet):
    points=[[Q(v) for v in row] for row in packet['vertices']]
    total=Q(0)
    for a,b,c in packet['faces']:
        p,q,r=points[a],points[b],points[c]
        total+=sum(p[i]*(q[(i+1)%3]*r[(i+2)%3]-q[(i+2)%3]*r[(i+1)%3]) for i in range(3))/6
    return total

class ExactUnionRecipe(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name); self.count=0
        self.executable=os.environ['EXACT_BOOLEAN_PROBE']

    def boolean(self,a,b,operation):
        self.count+=1; input_path=self.root/f'{self.count}.txt'; output=self.root/f'{self.count}.json'
        with input_path.open('x') as stream:
            for packet in (a,b):
                stream.write(f"{len(packet['vertices'])} {len(packet['faces'])}\n")
                for row in packet['vertices']:stream.write(' '.join(row)+'\n')
                for row in packet['faces']:stream.write(' '.join(map(str,row))+'\n')
        result=subprocess.run([self.executable,str(input_path),str(output),operation],
            text=True,capture_output=True,timeout=60)
        self.assertEqual(result.returncode,0,result.stderr)
        return json.loads(output.read_text())

    def test_union_of_overlapping_boxes_has_hand_computed_volume(self):
        a=box((-1,-1,-1),(1,1,1));b=box((0,-1,-1),(2,1,1))
        union=self.boolean(a,b,'union')
        self.assertEqual(exact_volume(union),12)
        self.assertFalse(self.boolean(union,box((-1,-1,-1),(2,1,1)),'difference')['faces'])
        self.assertFalse(self.boolean(box((-1,-1,-1),(2,1,1)),union,'difference')['faces'])

    def test_one_union_subtraction_equals_sequential_subtractions_and_keeps_shell(self):
        a=box((-1,-1,-1),(1,1,1));b=box((0,-1,-1),(2,1,1))
        inner=box((-2,-2,-2),(2,2,2));outer=box((-3,-3,-3),(3,3,3))
        union=self.boolean(a,b,'union')
        design=self.boolean(inner,union,'difference')
        sequential=self.boolean(self.boolean(inner,a,'difference'),b,'difference')
        self.assertEqual(exact_volume(design),52)
        self.assertFalse(self.boolean(design,sequential,'difference')['faces'])
        self.assertFalse(self.boolean(sequential,design,'difference')['faces'])
        assembled=self.boolean(outer,union,'union')
        fixed=self.boolean(assembled,design,'difference')
        self.assertEqual(exact_volume(assembled),216)
        self.assertEqual(exact_volume(fixed),164)
        self.assertFalse(self.boolean(fixed,design,'intersection')['faces'])

    def test_balanced_union_includes_every_part_and_matches_analytic_box(self):
        try:
            from balanced_exact_input import balanced_union
        except ImportError:
            balanced_union=None
        self.assertTrue(callable(balanced_union),'balanced exact construction is missing')
        parts=[box((i-1,-1,-1),(i+1,1,1)) for i in range(4)]
        union=balanced_union(parts,lambda a,b:self.boolean(a,b,'union'))
        self.assertEqual(exact_volume(union),20)
        reference=box((-1,-1,-1),(4,1,1))
        self.assertFalse(self.boolean(union,reference,'difference')['faces'])
        self.assertFalse(self.boolean(reference,union,'difference')['faces'])

if __name__=='__main__':unittest.main()
