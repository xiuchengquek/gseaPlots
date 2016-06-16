import argparse

class gseaParsers:
    @staticmethod
    def parse_rnk_file(filename):
        results = []
        with open(filename) as f:
            for line in f:
                line = line.strip()
                if not line.startswith('#'):
                    fields = line.split('\t')
                    fields[1] = float(fields[1])
                    results.append(tuple(fields))
        return results

    @staticmethod
    def parse_chip_file(filename):
        gene_id_to_probe_map = {}
        with open(filename) as f:
            next(f)
            for line in f:
                line = line.strip()
                fields = line.split('\t')
                gene_id_to_probe_map[fields[0]] = fields[1]
        return gene_id_to_probe_map

    @staticmethod
    def parse_single_gene_sets(filename):
        gene_set = {}
        with open(filename) as f:
            for line in f :
                line = line.strip()
                field = line.split('\t')
                gene_set[field[0]] = field[2:]
        return gene_set


class rankedList:
    def __init__(self, rnk_list, gene_chip):
        self.rnk_list = []
        self.gene_chip = gene_chip
        self.genes_set = {}
        self.enrichment_score = {}


        for x in rnk_list:
            genes = x[0]
            logFC = x[1]
            corrected_name = gene_chip.get(genes, genes)
            self.rnk_list.append((corrected_name, logFC) )

        self.rnk_list.sort(key= lambda x : x[1], reverse=True)

    def add_gene_set(self, genes_list, id ):
        self.genes_set[id] = genes_list

    def match_gene_set(self, id):
        self.enrichment_score[id] = []

        rolling_enrichment = 0
        sum_penalty = 0
        sum_hit = 0

        rnk_list = self.rnk_list
        genes_set = self.genes_set[id]

        penalty = self._calculate_gravity(genes_set , rnk_list)
        total_rank_score = self._calculate_sum_rank(genes_set, rnk_list)

        for rank, gene in enumerate(rnk_list):
            gene_id = gene[0]
            fold_change = gene[1]
            if gene_id not in genes_set:
                sum_penalty += penalty
            else :
                sum_hit += ( abs(fold_change) / total_rank_score )

            rolling_enrichment = sum_hit - sum_penalty
            self.enrichment_score[id].append((gene_id, rolling_enrichment))


    def write_enrichment_to_file(self, gene_set, file):
        if gene_set not in self.enrichment_score:
            pass
        else :
            with open(file, '\w+') as f :
                enrichment = self.enrichment_score[gene_set]
                for x in enrichment:
                    gene_id = x[0]
                    enrichment_value = x[1]
                    f.write("%s\t%s\n" % (gene_id, enrichment_value))

    def _calculate_sum_rank(self, genes_set, rnk_list):
        rnk_dict = {x[0] : x[1] for x in rnk_list}
        sum_rank = 0
        for x in genes_set:
            rank_score = rnk_dict.get(x, 0)
            sum_rank += abs(rank_score)
        return sum_rank

    def _calculate_gravity(self, genes_set, rnk_list):
        return 1 / (  len(rnk_list) - len(genes_set) )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate Cumulative Data for GSEA plot')
    parser.add_argument('-f', '--rnk' ,help='rnk file')
    parser.add_argument('-g', '--geneset' ,help='geneset gmt file')
    parser.add_argument('-c', '--chipfile', help='gene chip file')
    parser.add_argument('-o', '--output' , help='output file')

    args = parser.parse_args()

    rnk_list = gseaParsers.parse_rnk_file(args.rnk)
    gene_set = gseaParsers.parse_single_gene_sets(args.geneset)
    chip_file = gseaParsers.parse_chip_file(args.chipfile)

    output_file = args.output


    gene_set_id = list(gene_set.keys())[0]
    gene_set_values = gene_set[gene_set_id]


    rl = rankedList(rnk_list, chip_file)
    rl.add_gene_set(gene_set_values, gene_set_id)
    rl.match_gene_set(gene_set_id)
    rl.write_enrichment_to_file(gene_set_id, output_file)






























