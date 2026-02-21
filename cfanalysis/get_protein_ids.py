"""
Author: Patryk Jarnot (2026)
"""

import os


class ProteinIds:

    def __init__(self):
        self.gbsc_clusters_dir = None
        self.protein_id_path = None
        self.protein_acs = set()

    def run(self):
        if self.gbsc_clusters_dir is None:
            raise ValueError("gbsc dir not set")
        if self.protein_id_path is None:
            raise ValueError("output path not set")
        for file in os.listdir(self.gbsc_clusters_dir):
            with open(os.path.join(self.gbsc_clusters_dir, file), "r") as f:
                for line in f:
                    if not line.startswith(">"):
                        continue
                    items = line.split("|")
                    self.protein_acs.add(items[1])

        with open(self.protein_id_path, "w") as f:
            for uniprot_ac in self.protein_acs:
                f.write(f"{uniprot_ac}\n")
