import unittest
import SWAMP


class TestSequenceMasking(unittest.TestCase):


    def test_read_branchcodes(self):
        infile = 'example_dataset/branchcodes.txt'
        branch_codes = SWAMP.read_branchcodes(infile)

        # Branch codes come out like '5..7' -> (papio, colobus)
        branch_a = branch_codes['5..7']
        self.assertIn('papio', branch_a)
        self.assertIn('colobus', branch_a)

        # Single-species branch should have single-length array, e.g.
        # '6..2' -> (human)
        branch_b = branch_codes['6..2']
        self.assertEqual(len(branch_b), 1)


    def test_branch_error_check(self):
        seq_file = 'example_dataset/data/44/44.phy'
        seq_dict = SWAMP.read_phylip(seq_file)

        branch_file = 'example_dataset/branchcodes.txt'
        branch_codes = SWAMP.read_branchcodes(branch_file)

        # No error raised when valid branch codes provided.
        SWAMP.branch_error_check(branch_codes, seq_dict)

        # Error is raised when branchcodes is not present
        # (e.g. no file provided by user on command line)
        with self.assertRaises(ValueError) as cm:
            SWAMP.branch_error_check(None, seq_dict)

        # Mess up one of the branch_codes, adding an unknown species
        branch_codes['5..7'] += tuple(['gorilla'])
        # Error is raised when species isn't found in seq_dict.
        with self.assertRaises(ValueError) as cm:
            SWAMP.branch_error_check(branch_codes, seq_dict)


    def test_read_phylip(self):
        infile = 'example_dataset/data/101/101.phy'
        seq_dict = SWAMP.read_phylip(infile)

        # read_phylip returns a dict
        self.assertTrue(type(seq_dict) is dict)

        # Dict contains 4 seqs
        self.assertEqual(len(seq_dict.keys()), 4)

        # Sequence length check.
        self.assertEqual(len(seq_dict['pongo']), 1413)

        # Try another file.
        infile = 'example_dataset/data/44/44.phy'
        seq_dict = SWAMP.read_phylip(infile)

        self.assertEqual(len(seq_dict['pongo']), 1905)


    def test_write_phylip(self):
        infile = 'example_dataset/data/44/44.phy'
        seq_dict = SWAMP.read_phylip(infile)

        # Create some fake codons to mask.
        codons_to_mask = {}
        for i in range(2, 13):
            codons_to_mask[i] = ['pongo', 'homo']

        # Write file...
        masked_dict = SWAMP.mask_codons(seq_dict, codons_to_mask)
        SWAMP.print_masked_phyfile(infile, masked_dict)
        # Read it back
        masked_file = 'example_dataset/data/44/44_masked.phy'
        new_masked_dict = SWAMP.read_phylip(masked_file)

        pongo_seq = new_masked_dict['pongo']
        human_seq = new_masked_dict['homo']
        colobus_seq = new_masked_dict['colobus']

        # Make sure we see the masked seqs
        self.assertTrue(('NNN' * 10) in pongo_seq)
        self.assertTrue(('NNN' * 10) in human_seq)
        self.assertFalse(('NNN' * 10) in colobus_seq)


    def test_sliding_window_scan(self):
        infiles = ['example_dataset/data/44/44.phy']
        threshold = 1
        windowsize = 20
        interscan = False

        branch_file = 'example_dataset/branchcodes.txt'
        branch_codes = SWAMP.read_branchcodes(branch_file)

        # Run a sliding window scan on this single file.
        SWAMP.sliding_window_scan(infiles, threshold, windowsize,
                                  interscan, branch_codes)

        # Check that the masked file exists and contains 'NNN's
        masked_file = 'example_dataset/data/44/44_masked.phy'
        masked_dict = SWAMP.read_phylip(masked_file)
        pongo_seq = masked_dict['pongo']
        self.assertTrue('NNN' in pongo_seq)

        # Increase threshold, ensure NNNs are gone.
        threshold = 10
        SWAMP.sliding_window_scan(infiles, threshold, windowsize,
                                  interscan, branch_codes)
        masked_dict = SWAMP.read_phylip(masked_file)
        pongo_seq = masked_dict['pongo']
        self.assertFalse('NNN' in pongo_seq)


    def test_read_rst(self):
        # Extract branch information from rst file
        infile = 'example_dataset/data/44/44.phy'
        branches = SWAMP.read_rst(infile)

        self.assertEqual(len(branches['5..7']), 6)
        self.assertEqual(len(branches['5..6']), 0)


    def test_interscan(self):
        infile = 'example_dataset/data/44/44.phy'

        threshold = 1
        windowsize = 10
        interscan = False

        branch_file = 'example_dataset/branchcodes.txt'
        branch_codes = SWAMP.read_branchcodes(branch_file)

        # Mask without interscan... fewer masked codons
        result = SWAMP.sliding_window_scan_file(infile, threshold, windowsize,
                                                interscan, branch_codes)
        self.assertEqual(result['masked_column_count'], 215)

        # Mask with interscan... more masked codons
        interscan = True
        result = SWAMP.sliding_window_scan_file(infile, threshold, windowsize,
                                                interscan, branch_codes)
        self.assertEqual(result['masked_column_count'], 301)


if __name__ == '__main__':
    unittest.main()
