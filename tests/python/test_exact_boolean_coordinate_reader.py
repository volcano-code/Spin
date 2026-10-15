"""Coordinate binding before orientation must not certify material validity."""
import csv
from fractions import Fraction
from pathlib import Path
import tempfile
import unittest
import numpy as np
from make_it_stand.verification.rational_text import fraction_to_rational_text, fraction_from_rational_text

from make_it_stand.verification import faster_material as module

class CoordinateReaderTests(unittest.TestCase):
    def fixture(self,epsilon):
        exact=((Fraction(0),Fraction(0),Fraction(1)),
               (Fraction(1,2),Fraction(0),1+epsilon),
               (Fraction(1),Fraction(0),Fraction(1)),
               (Fraction(0),Fraction(1),Fraction(0)))
        points=np.asarray(exact,dtype=float)
        faces=np.asarray(((0,2,1),(0,1,3),(0,3,2),(1,2,3)))
        return exact,points,faces

    def write(self,path,exact,points,*,corrupt=None):
        headers=('vertex_index','status','exact_x','exact_y','exact_z',
                 'memory_x','memory_y','memory_z','obj_x','obj_y','obj_z')
        with path.open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=headers,delimiter='\t');writer.writeheader()
            for index,point in enumerate(exact):
                row={'vertex_index':str(index),'status':'MATCHED'}
                for axis,column in enumerate('xyz'):
                    row['exact_'+column]=fraction_to_rational_text(point[axis])
                    row['memory_'+column]=row['obj_'+column]=repr(float(points[index,axis]))
                if index==1 and corrupt:corrupt(row)
                writer.writerow(row)

    def test_bound_coordinates_are_available_before_region_validation(self):
        reader=getattr(module,'read_exact_boolean_sidecar_coordinates',None)
        self.assertTrue(callable(reader),'bound coordinate-only reader is missing')
        exact,points,faces=self.fixture(Fraction(1,2**60))
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'points.tsv';self.write(path,exact,points)
            rows,audit=reader(points,path)
            self.assertEqual(tuple(tuple(map(fraction_from_rational_text,row)) for row in rows),exact)
            self.assertNotIn('rational_geometry_identity',audit)
            self.assertNotIn('boundary_valid',audit)
            snapshot,complete=module.mesh_snapshot_from_exact_boolean_sidecar(points,faces,path,return_audit=True)
            self.assertTrue(snapshot.exact_region().facts['valid'])
            self.assertEqual(complete['sidecar_sha256'],audit['sidecar_sha256'])
            self.assertEqual(complete['rational_geometry_identity'],snapshot.identity)

    def test_bound_coordinates_never_upgrade_exact_collinearity(self):
        reader=getattr(module,'read_exact_boolean_sidecar_coordinates',None)
        self.assertTrue(callable(reader))
        exact,points,faces=self.fixture(Fraction(0))
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'points.tsv';self.write(path,exact,points)
            rows,_=reader(points,path)
            self.assertEqual(len(rows),4)
            with self.assertRaises(ValueError):
                module.mesh_snapshot_from_exact_boolean_sidecar(points,faces,path)

    def test_reader_preserves_sidecar_rejection_contract(self):
        reader=getattr(module,'read_exact_boolean_sidecar_coordinates',None)
        self.assertTrue(callable(reader))
        exact,points,_=self.fixture(Fraction(1,2**60))
        for key,value in (('vertex_index','0'),('status','AMBIGUOUS'),
                          ('memory_z','1.1'),('obj_z','1.1'),('exact_z','1001/1000')):
            with self.subTest(key=key),tempfile.TemporaryDirectory() as directory:
                path=Path(directory)/'points.tsv'
                self.write(path,exact,points,corrupt=lambda row:row.update({key:value}))
                with self.assertRaises(ValueError):reader(points,path)

if __name__=='__main__':unittest.main()

