import unittest
from analyze_failure_tape import source_report

class FailureTapeAnalysis(unittest.TestCase):
    def test_exact_tangency_is_not_negative_orientation_evidence(self):
        packet=dict(normalized_coefficients=[-1,0,0,0,1,1,1,0,0,0],snapshots=[dict(
            source_vertices=[[0,0,0],[0,1,0],[0,0,1]],source_faces=[[0,1,2]],
            source_face_is_cap=[0])])
        report=source_report(packet)
        self.assertEqual(report['negative_faces'],[])
        self.assertEqual(report['positive_faces'],[])
        self.assertEqual(report['zero_faces'],[0])

if __name__=='__main__':unittest.main()
