"""
Tests for the MSI sites reader.
"""
import unittest
from pathlib import Path
import tempfile

from pygnome.parsers.msi import MsiSitesReader, MsiSiteRecord


class TestMsiSitesReader(unittest.TestCase):
    """Test cases for the MSI sites reader."""
    
    def setUp(self):
        """Set up test data."""
        # Create a temporary MSI sites file
        self.temp_file = tempfile.NamedTemporaryFile(delete=False, mode='w')
        self.temp_file.write(
            "chromosome\tlocation\trepeat_unit_length\trepeat_unit_binary\trepeat_times\t"
            "left_flank_binary\tright_flank_binary\trepeat_unit_bases\tleft_flank_bases\tright_flank_bases\n"
            "chr1\t10485\t4\t149\t3\t150\t685\tGCCC\tAGCCG\tGGGTC\n"
            "chr1\t10629\t2\t9\t3\t258\t409\tGC\tCAAAG\tCGCGC\n"
            "chr1\t10652\t2\t2\t3\t665\t614\tAG\tGGCGC\tGCGCG\n"
            "chr2\t20000\t3\t100\t4\t200\t300\tCAT\tGGGGG\tAAAAA\n"
        )
        self.temp_file.close()
        self.file_path = Path(self.temp_file.name)
    
    def tearDown(self):
        """Clean up test data."""
        # Remove the temporary file
        Path(self.file_path).unlink()
    
    def test_reader_as_context_manager(self):
        """Test using the reader as a context manager."""
        records = []
        with MsiSitesReader(file_path=self.file_path) as reader:
            for record in reader:
                records.append(record)
        
        # Check that we got the expected number of records
        self.assertEqual(len(records), 4)
        
        # Check the first record
        self.assertEqual(records[0].chromosome, "chr1")
        self.assertEqual(records[0].location, 10485)
        self.assertEqual(records[0].repeat_unit_length, 4)
        self.assertEqual(records[0].repeat_unit_binary, 149)
        self.assertEqual(records[0].repeat_times, 3)
        self.assertEqual(records[0].left_flank_binary, 150)
        self.assertEqual(records[0].right_flank_binary, 685)
        self.assertEqual(records[0].repeat_unit_bases, "GCCC")
        self.assertEqual(records[0].left_flank_bases, "AGCCG")
        self.assertEqual(records[0].right_flank_bases, "GGGTC")
    
    def test_reader_as_iterator(self):
        """Test using the reader as an iterator."""
        reader = MsiSitesReader(self.file_path)
        try:
            records = list(reader)
            
            # Check that we got the expected number of records
            self.assertEqual(len(records), 4)
            
            # Check the second record
            self.assertEqual(records[1].chromosome, "chr1")
            self.assertEqual(records[1].location, 10629)
            self.assertEqual(records[1].repeat_unit_bases, "GC")
            self.assertEqual(records[1].repeat_times, 3)
        finally:
            reader.close()
    
    def test_fetch(self):
        """Test fetching records in a specific region."""
        with MsiSitesReader(self.file_path) as reader:
            # Fetch records in chr1:10600-10700
            records = list(reader.fetch("chr1", 10600, 10700))
            
            # Check that we got the expected number of records
            self.assertEqual(len(records), 2)
            
            # Check the records
            self.assertEqual(records[0].chromosome, "chr1")
            self.assertEqual(records[0].location, 10629)
            self.assertEqual(records[1].chromosome, "chr1")
            self.assertEqual(records[1].location, 10652)
    
    def test_record_properties(self):
        """Test the properties of MSI site records."""
        with MsiSitesReader(self.file_path) as reader:
            records = list(reader)
            
            # Check end_location
            self.assertEqual(records[0].end_location, 10485 + (4 * 3) - 1)  # location + (repeat_unit_length * repeat_times) - 1
            
            # Check sequence
            self.assertEqual(records[0].sequence, "GCCC" * 3)
            
            # Check full_sequence
            self.assertEqual(records[0].full_sequence, "AGCCG" + "GCCC" * 3 + "GGGTC")


if __name__ == "__main__":
    unittest.main()