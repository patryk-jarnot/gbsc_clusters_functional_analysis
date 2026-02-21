from cfanalysis.download_go import DownloadGo
from cfanalysis.get_protein_ids import ProteinIds
from cfanalysis.gbsc_functional_analysis import FunctionalAnalysis

from optparse import OptionParser

class CfPipeline:
    def __init__(self, gbsc_clusters_dir, protein_id_path):
        self.gbsc_clusters_dir = gbsc_clusters_dir
        self.protein_id_path = protein_id_path
        self.exclude_IEA = "no"
        self.aspect = "F"
        self.output_dir = "./gbsc_functional_results/"
        self.cache_file = "cache.sqlite"
        self.alpha = 0.05
        self.log_file = "gbsc_functional_analysis.log"

        self.protein_ids = ProteinIds()
        self.download_go = DownloadGo()
        self.functional_analysis = FunctionalAnalysis()

    def run(self):
        self.protein_ids.protein_id_path = self.protein_id_path
        self.protein_ids.gbsc_clusters_dir = self.gbsc_clusters_dir
        self.protein_ids.run()

        self.download_go.params.protein_id_path = self.protein_id_path
        self.download_go.params.exclude_IEA = self.exclude_IEA
        self.download_go.params.aspect = self.aspect
        self.download_go.params.output_dir = self.output_dir
        self.download_go.params.cache_file = self.cache_file
        self.download_go.run()

        self.functional_analysis.params.gbsc_clusters = self.gbsc_clusters_dir
        self.functional_analysis.params.alpha = self.alpha
        self.functional_analysis.params.output_dir = self.output_dir
        self.functional_analysis.params.log_file = self.log_file
        self.functional_analysis.run()


def main(cfpipeline):
    cfpipeline.run()


def get_options():
    parser = OptionParser(description="desc")
    parser.add_option("-p", "--protein_id_path", dest="protein_id_path", default=None,
                      help="List of proteins for annotations", metavar="FASTA")
    parser.add_option("-e", "--exclude-iea", dest="exclude_IEA", default="no",
                      help="Exclude GO terms with IEA? yes/no", metavar="STRING")
    parser.add_option("-s", "--aspect", dest="aspect", default="F",
                      help="Aspect of GO", metavar="STRING")
    parser.add_option('-o', '--output-dir', dest="output_dir", default='./gbsc_functional_results/',
                      help='Project directory')
    parser.add_option('-c', '--cache-file', dest="cache_file", default='cache.sqlite',
                      help='Cache directory for data downloaded')
    parser.add_option("-g", "--gbsc-clusters", dest="gbsc_clusters_dir", default=None,
                      help="Path to the directory with GBSC clusters", metavar="DIR")
    parser.add_option("-a", "--alpha", dest="alpha", default=0.05, type="float",
                      help="Threshold of test significance", metavar="FLOAT")
    parser.add_option('-l', '--log-file', dest="log_file", default='gbsc_functional_analysis.log',
                      help='Log file name')
    options, args = parser.parse_args()

    return options, args

"""

        self.gbsc_clusters_dir = gbsc_clusters_dir
        self.protein_id_path = protein_id_path
        self.exclude_IEA = "no"
        self.aspect = "F"
        self.output_dir = "./gbsc_functional_results/"
        self.cache_file = "cache.sqlite"
        self.alpha = 0.05
        self.log_file = "gbsc_functional_analysis.log"
        
"""

if __name__ == "__main__":
    options, args = get_options()
    cfpipeline = CfPipeline(options.gbsc_clusters_dir, options.protein_id_path)
    cfpipeline.exclude_IEA = options.exclude_IEA
    cfpipeline.aspect = options.aspect
    cfpipeline.output_dir = options.output_dir
    cfpipeline.cache_file = options.cache_file
    cfpipeline.alpha = options.alpha
    cfpipeline.log_file = options.log_file

    # try:
    try:
        main(cfpipeline)
    except KeyboardInterrupt:
        print("Shutdown requested...exiting")
        cfpipeline.download_go.cache.cursor.close()
        cfpipeline.download_go.cache.connection.close()
    except Exception as e:
        cfpipeline.download_go.cache.cursor.close()
        cfpipeline.download_go.cache.connection.close()
        raise e
