# Series of functions to modify protein sequences:
# 
# - leave
# - remove SP via SignalP 3
# - randomised n times
# - reversed
# - SP placed at end Seq, Cterm
# - SP randomised n times
# 
# 
test_protein = "MAPAVAGGGRVRSNEPVAAAAAAPAASGKPCVCGFQVCACTGSAAVASAASSLDMDIVAMGQIGAVNDESWVGVELGEDGETDESGAAVDDRPVFRTEKIKGVLLHPYRVLIFVRL"

import random

def protein_reverse(protein):
    
    return protein[::-1]
protein_reverse(test_protein)

def protein_random(protein):
    p = list(protein)
    random.shuffle(p)
    return ''.join(p)

def protein_sp_remove(protein, drange):
    """
    Return protein sequecne with SP removed. Given that SignalP will output a D-score and position, e.g.  1-29 
    then we can split on these and use for 
    """
    sp_start, sp_end = drange.split("-")
    return protein[int(sp_end):]    



def protein_sp_cterm(protein, drange):
    """
    Return protein sequence with SP placed at end
    """
    sp_start, sp_end = drange.split("-")
    
    return protein[int(sp_end):] + protein[:int(sp_end)]  

def protein_sp_random(protein, drange):
    """
    Return protein sequence with SP placed at end
    """
    sp_start, sp_end = drange.split("-")
    sp_only = protein[:int(sp_end)]
    sp = list(sp_only)
    random.shuffle(sp)
    
    return ''.join(sp) + protein[int(sp_end):]

import unittest

class ProteinVariationTest(unittest.TestCase):
    
    test_protein = "MAPAVAGGGRVRSNEPVAAAAAAPAASGKPCVCGFQVCACTGSAAVASAASSLDMDIVAMGQIGAVNDESWVGVELGEDGETDESGAAVDDRPVFRTEKIKGVLLHPYRVLIFVRL"
    test_range = "1-29"
    test_protein_reverse = "LRVFILVRYPHLLVGKIKETRFVPRDDVAAGSEDTEGDEGLEVGVWSEDNVAGIQGMAVIDMDLSSAASAVAASGTCACVQFGCVCPKGSAAPAAAAAAVPENSRVRGGGAVAPAM"
    test_protein_sp_remove = "PCVCGFQVCACTGSAAVASAASSLDMDIVAMGQIGAVNDESWVGVELGEDGETDESGAAVDDRPVFRTEKIKGVLLHPYRVLIFVRL"
    test_protein_sp = "MAPAVAGGGRVRSNEPVAAAAAAPAASGK"
    test_protein_sp_cterm  = "PCVCGFQVCACTGSAAVASAASSLDMDIVAMGQIGAVNDESWVGVELGEDGETDESGAAVDDRPVFRTEKIKGVLLHPYRVLIFVRLMAPAVAGGGRVRSNEPVAAAAAAPAASGK" 
    
    def test_protein_reverse(self):
        test_protein = "MAPAVAGGGRVRSNEPVAAAAAAPAASGKPCVCGFQVCACTGSAAVASAASSLDMDIVAMGQIGAVNDESWVGVELGEDGETDESGAAVDDRPVFRTEKIKGVLLHPYRVLIFVRL"
        test_protein_reverse = "LRVFILVRYPHLLVGKIKETRFVPRDDVAAGSEDTEGDEGLEVGVWSEDNVAGIQGMAVIDMDLSSAASAVAASGTCACVQFGCVCPKGSAAPAAAAAAVPENSRVRGGGAVAPAM"
        self.assertEqual(test_protein_reverse, protein_reverse(test_protein))
        
    def test_protein_random(self):
        test_protein = "MAPAVAGGGRVRSNEPVAAAAAAPAASGKPCVCGFQVCACTGSAAVASAASSLDMDIVAMGQIGAVNDESWVGVELGEDGETDESGAAVDDRPVFRTEKIKGVLLHPYRVLIFVRL"

        # Minimal check that legnth is correct, and not equal to original
        self.assertEqual(len(test_protein), len(protein_random(test_protein)))
        self.assertNotEqual(test_protein, protein_random(test_protein))

    def test_protein_sp_remove(self):
        test_protein = "MAPAVAGGGRVRSNEPVAAAAAAPAASGKPCVCGFQVCACTGSAAVASAASSLDMDIVAMGQIGAVNDESWVGVELGEDGETDESGAAVDDRPVFRTEKIKGVLLHPYRVLIFVRL"
        test_range = "1-29"
        test_protein_sp_remove = "PCVCGFQVCACTGSAAVASAASSLDMDIVAMGQIGAVNDESWVGVELGEDGETDESGAAVDDRPVFRTEKIKGVLLHPYRVLIFVRL"
        
        self.assertEqual(test_protein_sp_remove, protein_sp_remove(test_protein,test_range))
                         
    def test_protein_sp_cterm(self):
        test_protein = "MAPAVAGGGRVRSNEPVAAAAAAPAASGKPCVCGFQVCACTGSAAVASAASSLDMDIVAMGQIGAVNDESWVGVELGEDGETDESGAAVDDRPVFRTEKIKGVLLHPYRVLIFVRL"
        test_range = "1-29"
        test_protein_sp_cterm  = "PCVCGFQVCACTGSAAVASAASSLDMDIVAMGQIGAVNDESWVGVELGEDGETDESGAAVDDRPVFRTEKIKGVLLHPYRVLIFVRLMAPAVAGGGRVRSNEPVAAAAAAPAASGK" 
    
        self.assertEqual(test_protein_sp_cterm, protein_sp_cterm(test_protein,test_range))
                         
    def test_protein_sp_random(self):
        test_protein = "MAPAVAGGGRVRSNEPVAAAAAAPAASGKPCVCGFQVCACTGSAAVASAASSLDMDIVAMGQIGAVNDESWVGVELGEDGETDESGAAVDDRPVFRTEKIKGVLLHPYRVLIFVRL"
        test_range = "1-29"

        # Minimal check that legnth is correct, and not equal to original
        self.assertEqual(len(test_protein), len(protein_sp_random(test_protein,test_range)))    
        self.assertNotEqual(test_protein, protein_sp_random(test_protein,test_range))

                         

suite = unittest.TestLoader().loadTestsFromTestCase(ProteinVariationTest)
unittest.TextTestRunner(verbosity=2).run(suite)
