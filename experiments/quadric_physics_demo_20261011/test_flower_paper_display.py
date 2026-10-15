import unittest


class DisplayUnitsTests(unittest.TestCase):
    def test_circle_uses_flower_request_units_instead_of_heart_units(self):
        from flower_paper_display import circle_in_mm
        center,radius=circle_in_mm({'circle_center':[.128,-.373],'circle_radius':.394},
                                   {'meters_per_model_unit':.02})
        self.assertAlmostEqual(center[0],2.56)
        self.assertAlmostEqual(center[1],-7.46)
        self.assertAlmostEqual(radius,7.88)

    def test_missing_or_invalid_units_cannot_be_silently_assumed(self):
        from flower_paper_display import circle_in_mm
        for value in (0.,float('nan'),-.02):
            with self.assertRaises(ValueError):
                circle_in_mm({'circle_center':[0.,0.],'circle_radius':1.},
                             {'meters_per_model_unit':value})


if __name__=='__main__':unittest.main()
