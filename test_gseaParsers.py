import unittest
from unittest import mock
from enrichmentScore import gseaParsers, rankedList



class TestGseaParser(unittest.TestCase):

    def setUp(self):
        self.rnk_file_content = [
            '#',
            'geneB\t2',
            'geneA\t3',
            'geneC\t-1'
        ]

        self.gene_chip_content = [
            "Probe Set ID\tGene Symbol\tGene Title\tAliases",
            "A1BG\tA1BG\talpha-1-B glycoprotein\tGABDKFZP686F0970ABGHYST2477A1B",
            "A2M\tA2M\talpha-2-macroglobulin\tDKFZP779B086S863-7ALPHA 2MFWP007CPAMD5",
            "A2ML1\tA2ML1\talpha-2-macroglobulin-like 1\tFLJ39129DKFZP686O1010F"
        ]

        self.gene_set_gmt = [
            'GenesetA\tDescriptionHere\tGeneA\tGeneB\tGeneC'
        ]

    def test_gseaParser_parse_rnk_file(self):
        m = unittest.mock.mock_open(read_data='\n'.join(self.rnk_file_content))
        m.return_value.__iter__ = lambda self: self
        m.return_value.__next__ = lambda self: self.readline()
        with mock.patch('builtins.open',m ):
            results = gseaParsers.parse_rnk_file('test.rnk')
            self.assertEqual(results, [('geneB', 2), ('geneA', 3), ('geneC', -1)])

    def test_parse_chip_file(self):
        m = unittest.mock.mock_open(read_data='\n'.join(self.gene_chip_content))
        m.return_value.__iter__ = lambda self: self
        m.return_value.__next__ = lambda self: self.readline()


        with mock.patch('builtins.open',m ):
            gene_id_to_probe_map = gseaParsers.parse_chip_file('test.chip')
            self.assertDictEqual(gene_id_to_probe_map,
                                 {   'A1BG' : 'A1BG', 'A2M' : 'A2M', 'A2ML1' : 'A2ML1' })

    def test_parse_single_gene_sets(self):
        m = unittest.mock.mock_open(read_data='\n'.join(self.gene_set_gmt))
        m.return_value.__iter__ = lambda self: self
        m.return_value.__next__ = lambda self: self.readline()
        with mock.patch('builtins.open',m ):
            gene_sets = gseaParsers.parse_single_gene_sets('test.chip')
            self.assertDictEqual(gene_sets , {'GenesetA' : ['GeneA','GeneB','GeneC']})


class TestRankedList(unittest.TestCase):
    def setUp(self):
        self.rnk_list = [('geneB', 2), ('geneA', 3), ('geneC', -1)]
        self.gene_set = {'GenesetA' : ['GeneA']}
        self.gene_id_to_probe_map = {   'geneA' : 'GeneA', 'geneB' : 'GeneB', 'geneC' : 'GeneC' }
        self.rank_obj = rankedList(self.rnk_list, self.gene_id_to_probe_map)
        self.rank_obj.add_gene_set( self.gene_set['GenesetA'],'GenesetA')


    def test_rankedList_initialization(self):
        ## notice that the rnklist and the non rank list are different

        self.assertListEqual(self.rank_obj.rnk_list, [('GeneA', 3),('GeneB', 2),('GeneC' , -1)])

    def test_rankedList_add_gene_set(self):
        self.assertDictEqual(
            self.rank_obj.genes_set, self.gene_set
        )

    def test_match_gene_set(self):
        gravity = self.rank_obj._calculate_gravity(self.gene_set, self.rank_obj.rnk_list)
        self.assertEqual(gravity, 0.5)

    def test_calculate_sum_rank(self):
        sum_rank = self.rank_obj._calculate_sum_rank(self.rank_obj.genes_set['GenesetA'] ,self.rank_obj.rnk_list)
        self.assertEqual(sum_rank, 3)

    def test_match_gene_set(self):
        self.rank_obj.match_gene_set('GenesetA')
        self.assertListEqual([1 , 0.5, 0] , self.rank_obj.enrichment_score['GenesetA'])

class TestRankedListNegative(unittest.TestCase):
    def setUp(self):
        self.rnk_list = [('geneB', -2), ('geneA', -3), ('geneC', 1)]
        self.gene_set = {'GenesetA' : ['GeneA']}
        self.gene_id_to_probe_map = {   'geneA' : 'GeneA', 'geneB' : 'GeneB', 'geneC' : 'GeneC' }
        self.rank_obj = rankedList(self.rnk_list, self.gene_id_to_probe_map)
        self.rank_obj.add_gene_set( self.gene_set['GenesetA'],'GenesetA')

    def test_match_gene_set(self):
        self.rank_obj.match_gene_set('GenesetA')
        self.assertListEqual([-0.5 , -1.0, 0.0] , self.rank_obj.enrichment_score['GenesetA'])




















































if __name__ == '__main__':
    unittest.main()
